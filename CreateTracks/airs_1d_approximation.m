clearvars
addpath('../common');
CoreVars = sampling_core_variables;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%extract lat,lon,z of each AIRS point, to allow model sampling
%
%properties are modified to simulate a 1D approximation - no rotation, minimal
%width
%
%store in daily files of a common format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Settings.Instrument = 'airs1d';
Settings.InDir      = CoreVars.Airs.Path;
Settings.OutDir     = [CoreVars.MasterPath,'/tracks/AIRS_1D/'];
Settings.PrsLevels  = CoreVars.Airs.HeightRange;
Settings.LatRange   = [-1,1].*90; 
Settings.TimeRange  = [datenum(200511,1,1),datenum(2011,12,31)];


% % % %get a mean density profile, to scale the weighting functions by
% % % Density = load('../common/saber_density_new.mat');
% % % Density.D = nanmean(Density.Results.Density,1);
% % % Density.Z = Density.Results.HeightScale;
% % % Density = rmfield(Density,{'Results'});


OldFile = '';
for iDay=Settings.TimeRange(1):1:Settings.TimeRange(2);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %generate file name
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  DayFile = [Settings.OutDir,'track_',Settings.Instrument,'_',num2str(iDay),'.mat'];
  
  
  if exist(DayFile); disp([DayFile,' already done']); continue; end  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %import data
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  [y,m,d] = datevec(iDay);
  FileName = ['AIRS.',sprintf('%04d',y),'.',sprintf('%02d',m),'.',sprintf('%02d',d),'.all.L1B.mat'];
  FilePath = [Settings.InDir,'/',sprintf('%04d',y),'/',FileName];
  clear y m d FileName
  if ~exist(FilePath); 
    clear FilePath; 
    disp([datestr(iDay),' input file not located'])
    continue; end %no data
  load(FilePath);
  clear FilePath;
  
  %discard bad data
  Bad = find(Airs.Latitude  == -9999 ...
           | Airs.Longitude == -9999 ...
           | Airs.Time      == -9999);
  Airs.Latitude( Bad) = NaN;
  Airs.Longitude(Bad) = NaN;
  Airs.Time(     Bad) = NaN;

  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %get the data we want
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

  %discard anything that doesn't overlap our lat range
  Granules = ones(size(Airs.GranulesStored));
  for i=1:1:numel(Granules);
    Lat = squeeze(Airs.Latitude(i,:,:));
    TheMax = max(Lat(:));
    TheMin = min(Lat(:));
    if TheMax < min(Settings.LatRange); Granules(i) = 0; end
    if TheMin > max(Settings.LatRange); Granules(i) = 0; end
    clear Lat TheMax TheMin
  end
  Granules = find(Granules == 1);
  Airs.Latitude  = Airs.Latitude( Granules,:,:);
  Airs.Longitude = Airs.Longitude(Granules,:,:);
  Airs.Time      = Airs.Time(     Granules,:,:);
  clear Granules
  
  %convert timestamps
  Airs.Time = cjw_time_airs2matlab(Airs.Time);
  
  %ok. duplicate out the data to have the right number of heights
  Latitude  = repmat(Airs.Latitude, 1,1,1,numel(Settings.PrsLevels));
  Longitude = repmat(Airs.Longitude,1,1,1,numel(Settings.PrsLevels));
  AirsTime  = repmat(Airs.Time,     1,1,1,numel(Settings.PrsLevels));
  Pressure  = permute(repmat(Settings.PrsLevels',1,size(Latitude,1),size(Latitude,2),size(Latitude,3)),[2 3 4 1]);
  

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %prepare the necessary information to reconstruct the granules from 1D data
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

  %prepare the necessary information to reconstruct the granules from 1D data
  g = 1:1:size(Latitude,1); %granule
  x = 1:1:size(Latitude,2); %cross-track
  y = 1:1:size(Latitude,3); %along-track
  z = 1:1:size(Latitude,4); %pressure levels
  [g,x,y,z] = ndgrid(g,x,y,z); %grid
  Recon.g = g(:);
  Recon.x = x(:);
  Recon.y = y(:);
  Recon.z = z(:);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %define the approximate weighting kernel
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
 
  %we're using L1 AIRS data, not L2.
  %so, we want to use WEIGHTING functions rather than KERNEL functions
  %these are not regular, but rather slanty
  %so load them up and use them
  
  %storage method:
  %keep each *unique* function
  %then map to them with NEGATIVE numbers in the weight array
  %the processing routine will recognise them and translate accordingly
 
  %first, store in one array the weights for all the levels we use
  if ~exist('Weight'); %i.e. first loop
    
    %get the pressure weighting
    ZFuncs = gong_channel_weights(Settings.PrsLevels);

        
    %map which goes with each MEASUREMENT pressure
    ZFuncs.Channels = Settings.PrsLevels;
    
    %also create an array which will index to this
    Weight.Z = NaN(numel(Pressure),1); %make sure the order stays the same! 
  end
  
  %fill in the indexing points
  for iPoint=1:1:numel(Pressure)
    [~,idx] = min(abs(ZFuncs.Channels - Pressure(iPoint)));
    Weight.Z(iPoint) = -idx;  
  end
  Weight.ZFuncs = ZFuncs;

  
  %the  horizontal weighting functions are approximately...
  Weight.X = ones(size(Recon.x)).*0.1./4; %along-LOS - 100m
  Weight.Y = ones(size(Recon.x)).*0.1./4; %across-LOS - 100m
  %(Hoffmann et al, AMT 2014)
  %see saber and hirdls routines for explanation of /4     
  
  
  %the vertical angle varies by distance off-axis geometrically
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %angle is defined as a/c/w from nadir
  
  ViewAngleZ = zeros(size(Recon.x)); %rows off-centre
  ViewAngleH = ViewAngleZ;
  

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %save!
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
  
  %put data into a line
  Prs  = Pressure(:);
  Time = AirsTime(:);
  Lat  = Latitude(:);
  Lon  = Longitude(:);  
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
