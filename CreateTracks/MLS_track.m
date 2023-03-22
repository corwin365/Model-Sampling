clearvars
addpath(genpath('../'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%produce OIF data for MLS
%
% Phoebe Noble 2023/02/09 - adapted from MLS script (Corwin Wright,
% c.wright@bath.ac.uk, 2023/02/03)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%dataset identifier
Settings.Instrument = 'MLS';

%where do the input files live?
Settings.InDir = [LocalDataDir,'MLS/'];

%geolocation - which data should we include?
%for all except HeightRange, we include any wholegranule including these
%for HeightRange, we will trim the granules in height to just this range
Settings.LatRange    = [-90,90];
Settings.LonRange    = [-180,180];
Settings.TimeRange   = [1,1].*datenum(2020,1,24);
Settings.HeightRange = [0,80]; %km

%path handling internal to routine
[~,CoreSettings] = sampling_core_v2(' ',' ',0,'GetSettings',true);
Settings.OutDir  = [CoreSettings.MasterPath,'/tracks/',Settings.Instrument,'/'];
clear CoreSettings



for iDay=Settings.TimeRange(1):1:Settings.TimeRange(2);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %generate file name
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  %create a directory to put the data in, if it doesn't exist
  if exist(Settings.OutDir,'dir') ~= 7; mkdir(Settings.OutDir); end

  %generate filename
  DayFile = [Settings.OutDir,'track_',Settings.Instrument,'_',num2str(iDay),'.mat'];
  
  
  if exist(DayFile); continue; end % if we've already made the track, don't do it again.
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %import data
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  %identify the relevant file and load it
  MLSData = extract_MLS_data(iDay,Settings.InDir);
  if MLSData.Error ~= 0; 
    %problem getting data - skip
    disp('Problem getting data')
    clear MLSData
    continue
  end
  

  
% % we don't need to calculate travel angle here unlike HIRDLS as we have
% it directly from the retrieval.
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   %travel angle!
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% 
%   %find the travel direction at each point
%   %assume one height level - should be pretty close
%   LatScale = MLSData.Lat(:,20);
%   LonScale = MLSData.Lon(:,20);
%   
%   %get the next point
% %   NextLat = circshift(LatScale,[0,-1]);
% %   NextLon = circshift(LonScale,[0,-1]);
% 
%   NextLat = circshift(LatScale,-1);
%   NextLon = circshift(LonScale,-1);
% 
%   %this introduce one bad point - kill it
%   LatScale = LatScale(2:end); LonScale = LonScale(2:end);
%   NextLat  = NextLat( 2:end); NextLon  = NextLon( 2:end);
%   
%   %then find the travel direction at each point
%   Azimuth = azimuth(LatScale,LonScale,NextLat,NextLon)'; 
%   Azimuth = [NaN,Azimuth];
%   Azimuth(1) = Azimuth(2); %should be very close


  % MLS instrument looks forward on the AURA satellite therefore
  ViewAngleH = MLSData.LineOfSightAngle;

  %z viewing angle is zero - instrument is viewing from the limb
  ViewAngleZ = zeros(size(ViewAngleH));  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %get the data we want
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
  %extract prs scale 
  Prs = MLSData.Prs(1,:);
  
  %reduce down to just the pressure range we want
  InPrsRange = find(Prs >= h2p(max(Settings.HeightRange)) ...
                  & Prs <= h2p(min(Settings.HeightRange)));
  Prs = Prs(  InPrsRange);
  Lat = MLSData.Lat(:,InPrsRange);
  Lon = MLSData.Lon(:,InPrsRange);
  ViewAngleH = ViewAngleH(:,InPrsRange);
  ViewAngleZ = ViewAngleZ(:,InPrsRange);  
  
  clear InPrsRange
  
  %generate 30km lat array, for trimming later
  [~,idx] = min(abs(p2h(30)-Prs));
  LatScale = Lat(:,idx);
  clear idx
  
  %merge into point x height grid
  [Prs,Time] = meshgrid(Prs,MLSData.Time(:,1)); %we assumed time constant across a profile   earlier

% % %   %and trim data down to region
% % %   InLatRange = find(LatScale >= min(Settings.LatRange) ...
% % %                   & LatScale <= max(Settings.LatRange));
% % %   Lat  = Lat( InLatRange,:);              
% % %   Lon  = Lon( InLatRange,:);
% % %   Time = Time(InLatRange,:);
% % %   Prs  = Prs( InLatRange,:);
% % %   ViewAngle = ViewAngle(InLatRange,:);  
% % %   clear InLatRange LatScale

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %prepare the necessary information to reconstruct the granules from 1D data
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

  x = 1:1:size(Lat,1);  
  z = 1:1:size(Lat,2);
  [z,x] = meshgrid(z,x);
  Recon.x = single(x(:));
  Recon.z = single(z(:));
  clear x z

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %define the approximate weighting blob
  %numbers I'm using are ~2 stdev, so divide appropriately
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
  
  %approximate shape of blob
  ShapeX  = 165; %165km along-LOS
  ShapeY1  = 6;  %6km across-LOS below m'sphere
  ShapeY2  = 12;  %12km across-LOS below m'sphere
  
  %%% Calculate weightings based on resolutions:
  % X
  Weight.X = single(ones(size(Recon.x)).*ShapeX./4); %/2 for stdevs, additional /2 for width rather than distance from centre 
  
  % Y 
  %separate Y in t'sphere and above t'sphere (defined by pressure level 0.01hPa)
  Weight.Y = single(ones(size(Recon.x))); 
  LowY = find(Prs(:) > 0.01);
  Weight.Y(LowY) = ShapeY1./4;
  
  HighY = find(Prs(:) <= 0.01);
  Weight.Y(HighY) = ShapeY2./4;
  
  % Z
  % Z is defined by table 3.22.1 in MLS v5-0_data_quality_document
  % which gives the following pressure values:
  given_pressures = [0.00022, 0.00046, 0.001, 0.01, 0.1, 0.316, 1, 3.6, 10,...
    14.7, 31.6, 56.2, 100, 215, 261];
  % and their corresponding resolutions:
  given_horz_resolution = [12, 13,12 11, 6.4, 8.1, 6.8, 5.5, 4.1, 3.9, 3.6, 3.7, 4.6, ...
    3.5, 3.8];
    
  % convert to height
  given_heights = p2h(given_pressures);
  % interpolate to nearest neighbor for our pressures.
  ShapeZ = interp1(given_heights, given_horz_resolution, p2h(Prs(:)), 'nearest');
  
  Weight.Z = ShapeZ./4;
  

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %save!
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
   
  
  %put data into a line
  Prs  = Prs(:);
  Time = Time(:);
  Lat  = Lat(:);
  Lon  = Lon(:);  
  ViewAngleH = ViewAngleH(:);
  ViewAngleZ = ViewAngleZ(:);  
  
  %then into an array
  Track.Lat       = single(Lat);
  Track.Lon       = single(Lon);
  Track.Prs       = single(Prs);
  Track.Time      = Time; %needs to be double
  Track.ViewAngleH = single(ViewAngleH);
  Track.ViewAngleZ = single(ViewAngleZ);
  
  clear Lat Lon Prs Time
  
  
  %and save it
  save(DayFile,'Track','Recon','Weight');
  
  
  %tidy up, then done
  %clear Track DayFile
  disp([datestr(iDay),' complete'])
end
clear iDay
