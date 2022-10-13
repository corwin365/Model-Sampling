clearvars
addpath('../common');
CoreVars = sampling_core_variables;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%AIRS vertical resolution, but on a 2km x2km cross track grid of width 400km,
%at a 60deg angle to the vertical
%
%then store in granule-level files of a common format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Settings.Instrument = 'peterp';
Settings.InDir      = CoreVars.Airs3D.Path;
Settings.OutDir     = [CoreVars.MasterPath,'/tracks/PETERP/'];
Settings.PrsLevels  = CoreVars.Airs3D.HeightRange;
Settings.LatRange   = [0,90];
Settings.LonRange   = [0,180];
Settings.TimeRange  = datenum(2018,11,[1,14]);

%define the vertical resolution of the retrieval at each height
%taken from figure 5a of Hoffmann and Alexander 2009, cz=50km curve
%first column is altitude in km, second colum is resolution in km
Res = [18 ,10.25;...
       21 ,6.5;...
       24 ,6.75;...
       27 ,7.25;...
       30 ,7.5;...
       33 ,7.75;...
       36 ,8.75;...
       39 ,9;...
       42 ,9;...
       45 ,10;...
       48 ,10;...
       51 ,11;...
       54 ,13.5;...
       57 ,14.5;...
       60 ,13.25;...
       65 ,13.5;...
       70 ,19.5];
     
Res(:,2) = Res(:,2) ./ 2.355; %fwhm -> stdev
Res(:,2) = Res(:,2) ./ 2;     %half st-dev, as used in routine


OldFile = '';
for iDay=Settings.TimeRange(1):1:Settings.TimeRange(2);
  
  %find granules on this day in our geographic region
  List = find_airs_overpasses([1,1].*iDay,Settings.LonRange,Settings.LatRange,[20,0],1);
  List = List(:,2);
  
  %loop over them
  for jGranule=1:1:numel(List);
    iGranule = List(jGranule);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %generate file name
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    OutFile = [Settings.OutDir,'track_',Settings.Instrument,'_',num2str(iDay),'_g',sprintf('%03d',iGranule),'.mat'];
    
    if exist(OutFile) ~= 0;
      disp([datestr(iDay),', granule ',sprintf('%03d',iGranule),' already done'])
      continue
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %import data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %create lists of 1D, 2D and 3D fields, for later use
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    Fields.D1 = {};
    Fields.D2 = {'l1_time','l1_lon','l1_lat','ret_press'};
    Fields.D3 = {'ret_temp'};
    
    
    [y,~,~] = datevec(iDay);
    dayno = date2doy(iDay);
    
    DataDir = [Settings.InDir,sprintf('%04d',y),'/',sprintf('%03d',dayno)];
    File = wildcardsearch(DataDir,['*_',sprintf('%03d',iGranule),'.nc']);
    if numel(File) == 0;
      disp([datestr(iDay),', granule ',sprintf('%03d',iGranule),' NOT PRESENT'])
      continue
    end
    
    
    %load file
    Data = cjw_readnetCDF(File{1});
    Data.ret_press = h2p(Data.ret_z);

    
    %small bug in some of the input files
    %if Data.ret_z is all zeroes, set it to the default
    if sum(Data.ret_z(:)) == 0;
      Data.ret_z = [0;3;6;9;12;15;18;21;24;27;30;33;36;39;42;45;48;51;54;57;60;65;70;75;80;85;90];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %convert to Peter's suggested format
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %first, find the start and end of the middle row
    Start = [Data.l1_lat(45,  1),Data.l1_lon(45,  1)];
    End   = [Data.l1_lat(45,end),Data.l1_lon(45,end)];

    %hence, find the distance and azimuth separating them
    [dx,az] = distance(Start(1),Start(2),End(1),End(2));
    dx = deg2km(dx);

    %create a grid between these points
    Spacing = [2,2]; %along, across track
    Width   = [-1,1].*400; %km across-track
    x = 0:Spacing(1):dx;
    y = Width(1):Spacing(2):Width(2);
    [xi,yi] = meshgrid(x,y);

    %work out the distance and bearing of each point from the first element at centre-track
    r = quadadd(xi,yi);
    th = wrapTo180(atan2d(yi,xi)+az);

    %and hence work out the lat and lon of each point
    [Lat,Lon] = reckon(Start(1),Start(2),km2deg(r),th);
    clear Start End dz az Spacing Width x y xi yi th r

    %hence, replace the geolocation
    Data.l1_lat = Lat;
    Data.l1_lon = Lon;
    Data.l1_time = ones(size(Lon)).*mean(Data.l1_time,'all','omitnan'); %close enough!
    Data.ret_temp = NaN([numel(Data.ret_z),size(Data.l1_lat)]);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %discard unwanted heights
    InHeightRange = find(Data.ret_z >= p2h(CoreVars.Airs3D.HeightRange(1)) ...
                       & Data.ret_z <= p2h(CoreVars.Airs3D.HeightRange(2)));
    Data.ret_z    = Data.ret_z(InHeightRange);
    Data.ret_temp = Data.ret_temp(InHeightRange,:,:);
    clear InHeightRange
    
    %duplicate out fields that need higher dimensions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    sz = size(Data.ret_temp);
    Data.ret_z = repmat(Data.ret_z,1,sz(2),sz(3));
    
    for iField=1:1:numel(Fields.D1);
      Field = Data.(Fields.D1{iField});
      Field = permute(Field,[2,3,1]);
      Field = repmat(Field,sz(1),sz(2),1);
      Data.(Fields.D1{iField}) = Field;
    end
    
    for iField=1:1:numel(Fields.D2);
      Field = Data.(Fields.D2{iField});
      Field = permute(Field,[3,1,2]);
      Field = repmat(Field,sz(1),1,1);
      Data.(Fields.D2{iField}) = Field;
    end
    
    %all the fields are 3D now
    Fields.All = [Fields.D1,Fields.D2,Fields.D3,'ret_z'];
    
    %set bad values to NaN (default is -9999)
    for iField=1:1:numel(Fields.All);
      Field = Data.(Fields.All{iField});
      Field(Field == -9999) = NaN;
      Data.(Fields.All{iField}) = Field;
    end
    clear Field iField sz
    
    %also some negative times - get rid
    Data.l1_time(Data.l1_time < 0) = NaN;
    
    
    %convert timestamps
    Data.l1_time = datenum(2000,1,1,0,0,(Data.l1_time));
    
    %some files have messed-up timestamps where the date is wrong but the time is fine - fix these using the filename as a guide
    Data.l1_time = Data.l1_time-floor(Data.l1_time)+iDay;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %prepare the necessary information to reconstruct the granules from 1D data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %prepare the necessary information to reconstruct the time series from 1D data
    x = 1:1:size(Data.l1_lat,2); %cross-track
    y = 1:1:size(Data.l1_lat,3); %along-track
    z = 1:1:size(Data.l1_lat,1); %pressure levels
    [x,y,z] = ndgrid(x,y,z); %grid
    Recon.x = x(:);
    Recon.y = y(:);
    Recon.z = z(:);
    clear x y z
    
    %make the original data match in order
    for iField=1:1:numel(Fields.All);
      Field = Data.(Fields.All{iField});
      Field = permute(Field,[2,3,1]);
      Data.(Fields.All{iField}) = Field;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %define the approximate weighting kernel
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    %viewing angle in the horizontal is necessary, since it defines the angle of the z-splurging
    %we need to use a row-based approach, otherwise we'll get a jump at the end of each height level
    LatScale = Data.l1_lat; NextLat = circshift(LatScale,[0,-1,0]);
    LonScale = Data.l1_lon; NextLon = circshift(LonScale,[0,-1,0]);
    
    %this introduce one bad point - kill it
    LatScale(:,end,:) = NaN; LonScale(:,end,:) = NaN;
    NextLat( :,end,:) = NaN; NextLon( :,end,:) = NaN;
    
    %then find the travel direction at each point
    Azimuth = azimuth(LatScale,LonScale,NextLat,NextLon);
    
    %and replace the point we removed with the next value, which should be very
    %close to being correct
    Azimuth(:,end,:) = Azimuth(:,end-1,:);
    
    %finally, convert to a viewing angle
    ViewAngleH = wrapTo180(360-((360-Azimuth)))+45; %degrees c/w from N
    
    %smooth over minor jumps caused by finite precision on the geolocation
    ViewAngleH = smoothn(ViewAngleH,[5,5,1]);
    
    clear LatScale LonScale NextLat NextLon Azimuth

    %the  horizontal weighting functions are approximately...
    Weight.X = ones(size(Recon.x)).*13.5./(2.*2.355); %along-LOS
    Weight.Y = ones(size(Recon.x)).*13.5./(2.*2.355); %across-LOS
    %(Hoffmann et al, AMT 2014)
    
    %the vertical resolution is a function of altitude, as defined above
    Weight.Z = interp1(Res(:,1),Res(:,2),Data.ret_z(:),'linear','extrap');
    
% %     %the vertical angle varies by distance off-axis geometrically
% %     %angle is defined as a/c/w from nadir
% %     ViewAngleZ = mod(Recon.x,90) - (90./2); %rows off-centre
% %     ViewAngleZ = ViewAngleZ ./ max(abs(ViewAngleZ)); %normalised
% %     ViewAngleZ = ViewAngleZ .* 49.5; %degrees
    
    %for this study, the vertical angle is sublimb at 60deg from vertical
    ViewAngleZ = ones(size(Recon.x)).*60;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %save!
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %put data into a line
    Prs  = h2p(Data.ret_z(:));
    Time = Data.l1_time(:);
    Lat  = Data.l1_lat(:);
    Lon  = Data.l1_lon(:);
    ViewAngleH = ViewAngleH(:);
    ViewAngleZ = ViewAngleZ(:);
    
    %then into an array
    Track.Lat  = single(Lat);
    Track.Lon  = single(Lon);
    Track.Prs  = single(Prs);
    Track.Time = Time; %needs to be double
    Track.ViewAngleZ = single(ViewAngleZ);
    Track.ViewAngleH = single(ViewAngleH);
    clear Lat Lon Prs Time
    
    %and save it
    stop
    save(OutFile,'Track','Recon','Weight');
    
    
    
    %tidy up, then done
    clear Track DayFile
    disp([datestr(iDay),', granule ',sprintf('%03d',iGranule),' complete'])
  end; clear iGranule
end; clear iDay