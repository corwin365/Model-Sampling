clearvars
addpath(genpath('../'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%produce OIF data for 3D AIRS, as produced by the Hoffmann and Alexander 
% (2009) retrieval
%
%Corwin Wright, c.wright@bath.ac.uk, 2023/02/03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%dataset identifier
Settings.Instrument = 'AIRS3D';

%where do the input files live?
Settings.InDir = [LocalDataDir,'/AIRS/3d_airs'];

%geolocation - which data should we include?
%for all except HeightRange, we include any wholegranule including these
%for HeightRange, we will trim the granules in height to just this range
Settings.LatRange    = [-90,90];
Settings.LonRange    = [-180,180];
% Settings.TimeRange   = datenum(2018,11,5):1:datenum(2018,11,5)
Settings.TimeRange   = datenum(2020,1,20):1:datenum(2020,3,1);
Settings.HeightRange = [20,60]; %km

%path handling internal to routine
[~,CoreSettings] = sampling_core_v2(' ',' ',0,'GetSettings',true);
Settings.OutDir  = [CoreSettings.MasterPath,'/tracks/',Settings.Instrument,'/'];
clear CoreSettings

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameters used internally
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create lists of 1D, 2D and 3D fields, for setup use
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Fields.D1 = {};
Fields.D2 = {'l1_time','l1_lon','l1_lat','ret_press'};
Fields.D3 = {'ret_temp'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%main extraction-and-creation loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iDay=min(Settings.TimeRange):1:max(Settings.TimeRange);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  %find granules on this day in our geographic region
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  List = find_airs_overpasses(iDay+[0,1],Settings.LonRange,Settings.LatRange,[0,0],1);
  List = List(:,2);

  %loop over them
  for jGranule=1:1:240%numel(List);
    try
  % for jGranule=1:1:numel(List);
    iGranule = List(jGranule);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %generate output file name
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %create a directory to put the data in, if it doesn't exist
    if exist(Settings.OutDir,'dir') ~= 7; mkdir(Settings.OutDir); end

    %and generate the file name
    OutFile = [Settings.OutDir,'track_',Settings.Instrument,'_',num2str(iDay),'_g',sprintf('%03d',iGranule),'.mat'];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %import data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %find file holding the data for this day, and import it
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [y,~,~] = datevec(iDay);
    dayno = date2doy(iDay);
    
    DataDir = [Settings.InDir,'/',sprintf('%04d',y),'/',sprintf('%03d',dayno)];
    File = wildcardsearch(DataDir,['*_',sprintf('%03d',iGranule),'.nc']);
    if numel(File) == 0;
      disp([datestr(iDay),', granule ',sprintf('%03d',iGranule),' NOT PRESENT'])
      continue
    end
    
    %load file
    Data = rCDF(File{1});
    
    
    %tidy data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %small bug in some of the input files
    %if Data.ret_z is all zeroes, set it to the default
    if sum(Data.ret_z(:)) == 0;
      Data.ret_z = [0;3;6;9;12;15;18;21;24;27;30;33;36;39;42;45;48;51;54;57;60;65;70;75;80;85;90];
    end

    Data.ret_press = h2p(Data.ret_z);    

    %discard unwanted heights
    InHeightRange = find(Data.ret_z >= min(Settings.HeightRange) ...
                       & Data.ret_z <= max(Settings.HeightRange) );
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
    
    %convert timestamps to Matlab
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
    Weight.X = ones(size(Recon.x)).*13.5./(2.*2.355); %along-LOS. the factor of 2.355converts from FWHM to stdev, the factor of 2 from stdev to the stdev/2 used internally
    Weight.Y = ones(size(Recon.x)).*13.5./(2.*2.355); %across-LOS
    %(Hoffmann et al, AMT 2014)
    
    %compute vertical resolution, based on mean geographic location of granule
    [latmean,~] = meanm(Data.l1_lat(:),Data.l1_lon(:));
    RNight = airs_resolution(1,dayno,latmean,squeeze(Data.ret_z(1,1,:)))./(2.*2.355);
    RDay   = airs_resolution(0,dayno,latmean,squeeze(Data.ret_z(1,1,:)))./(2.*2.355);
    
    IsDay =  which_airs_retrieval(Data.l1_lon(:),Data.l1_lat(:),Data.l1_time(:),1);
    Weight.Z = Weight.X.*NaN;
    Weight.Z(IsDay == 1) = abs(interp1(squeeze(Data.ret_z(1,1,:)),  RDay,Data.ret_z(IsDay == 1))); 
    Weight.Z(IsDay == 0) = abs(interp1(squeeze(Data.ret_z(1,1,:)),RNight,Data.ret_z(IsDay == 0)));
    clear latmean RNight RDay IsDay

    %the vertical angle varies by distance off-axis geometrically
    %angle is defined as a/c/w from nadir
    ViewAngleZ = mod(Recon.x,90) - (90./2); %rows off-centre
    ViewAngleZ = ViewAngleZ ./ max(abs(ViewAngleZ)); %normalised
    ViewAngleZ = ViewAngleZ .* 49.5; %degrees

    %there is a small bug I need to pin down in the main code that has trouble dealing with ViewAngleZ *exactly* equal to 0
    %as a temporary patch that will have no meaningful effect on the results, set these valeus to a nonzero value
    ViewAngleZ(ViewAngleZ == 0) = 0.001;
    
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

    save(OutFile,'Track','Recon','Weight');
    
    
    
    %tidy up, then done
    clear Track DayFile
    disp([datestr(iDay),', granule ',sprintf('%03d',iGranule),' complete'])
    catch; end

  end; clear iGranule
end; clear iDay
