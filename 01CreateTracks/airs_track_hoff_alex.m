clearvars
addpath('../common');
CoreVars = sampling_core_variables;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%extract lat,lon,z of each AIRS point, to allow model sampling
%store in daily files of a common format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Settings.Instrument = 'airs3d';
Settings.InDir      = CoreVars.Airs3D.Path;
Settings.OutDir     = [CoreVars.MasterPath,'/tracks/AIRS_3D/'];
Settings.PrsLevels  = CoreVars.Airs3D.HeightRange;
Settings.LatRange   = [-1,1].*90; 
Settings.TimeRange  = CoreVars.Airs3D.TimeRange;


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

     
OldFile = '';
for iDay=Settings.TimeRange(1):1:Settings.TimeRange(2);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %generate file name
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  DayFile = [Settings.OutDir,'track_',Settings.Instrument,'_',num2str(iDay),'.mat'];
  
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
  Files = wildcardsearch(DataDir,'*.nc');
  NFiles = numel(Files);
  
  for iFile=1:1:NFiles
    
    %load file
    File = Files{iFile};
    Data = cjw_readnetCDF(File);
    
    %small bug in some of the input files
    %if Data.ret_z is all zeroes, set it to the default
    if sum(Data.ret_z(:)) == 0; 
      Data.ret_z = [0;3;6;9;12;15;18;21;24;27;30;33;36;39;42;45;48;51;54;57;60;65;70;75;80;85;90];
    end
    
    %remove unneeded fields
    Data = rmfield(Data,{'l1_nu','l1_rad','l2_time','l2_z','l2_lon','l2_lat','l2_press','l2_temp','l1_sat_z','l1_sat_lon','l1_sat_lat'});
    
    %discard unwanted heights
    InHeightRange = find(Data.ret_z >= p2h(CoreVars.Airs3D.HeightRange(1)) ...
                       & Data.ret_z <= p2h(CoreVars.Airs3D.HeightRange(2)));
    Data.ret_z    = Data.ret_z(InHeightRange);
    Data.ret_temp = Data.ret_temp(InHeightRange,:,:);
    
    %duplicate out the ones that need to be higher dimensions
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
    
    %also some negative times - get rid
    Data.l1_time(Data.l1_time < 0) = NaN;
    
    
    %store data
    if iFile==1;
      %create arrays
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
      AllData = Data;
    else  
      %concatenate the data to what we already have
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
      for iField=1:1:numel(Fields.All);
        Field = Fields.All{iField};
        AllData.(Field) = cat(3,AllData.(Field),Data.(Field));
      end
      
    end
    
    %done!
    clear Data sz iFile iField Field File InHeightRange y dayno

  end; clear iFile NFiles Files Fields

  %convert timestamps
  AllData.l1_time = cjw_time_airs2matlab(AllData.l1_time);
  

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %prepare the necessary information to reconstruct the granules from 1D data
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

  %prepare the necessary information to reconstruct the time series from 1D data
  x = 1:1:size(AllData.l1_lat,2); %cross-track
  y = 1:1:size(AllData.l1_lat,3); %along-track
  z = 1:1:size(AllData.l1_lat,1); %pressure levels
  [x,y,z] = ndgrid(x,y,z); %grid
  Recon.x = x(:);
  Recon.y = y(:);
  Recon.z = z(:);

  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %define the approximate weighting kernel
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     

  
  %viewing angle in the horizontal is zero, since we're a nadir sounder
  ViewAngleH = zeros(size(Recon.x));

  %the  horizontal weighting functions are approximately...
  Weight.X = ones(size(Recon.x)).*13.5./4; %along-LOS
  Weight.Y = ones(size(Recon.x)).*13.5./4; %across-LOS
  %(Hoffmann et al, AMT 2014)
  %see saber and hirdls routines for explanation of /4   
  
  %the vertical resolution is a function of altitude, as defined above
  Weight.Z = interp1(Res(:,1),Res(:,2),AllData.ret_z(:),'linear','extrap');
  
  %the vertical angle varies by distance off-axis geometrically
  %angle is defined as a/c/w from nadir
  ViewAngleZ = mod(Recon.y,135) - (135./2); %rows off-centre
  ViewAngleZ = ViewAngleZ ./ max(abs(ViewAngleZ)); %normalised
  ViewAngleZ = ViewAngleZ .* 49.5; %degrees


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %save!
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
  
  %put data into a line
  Prs  = h2p(AllData.ret_z(:));
  Time = AllData.l1_time(:);
  Lat  = AllData.l1_lat(:);
  Lon  = AllData.l1_lon(:); 
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
  save(DayFile,'Track','Recon','Weight');
  
  
  
  %tidy up, then done
  clear Track DayFile
  disp([datestr(iDay),' complete'])
end
clear iDay