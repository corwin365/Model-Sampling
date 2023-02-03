clearvars
addpath('../common');
CoreVars = sampling_core_variables;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%extract lat,lon,z of each HIRDLS point, to allow model sampling
%store in daily files of a common format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Settings.Instrument = 'hirdls';
Settings.InDir      = CoreVars.Hirdls.Path;
Settings.OutDir     = [CoreVars.MasterPath,'/tracks/HIRDLS/'];
Settings.PrsRange   = CoreVars.Hirdls.HeightRange; 
Settings.TimeRange  = [1,1].*datenum(2007,1,2);%CoreVars.Hirdls.TimeRange;

for iDay=Settings.TimeRange(1):1:Settings.TimeRange(2);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %generate file name
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  DayFile = [Settings.OutDir,'track_',Settings.Instrument,'_',num2str(iDay),'.mat'];
  
  
  if exist(DayFile); continue; end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %import data
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  %identify the relevant file and load it
  HirdlsData = extract_hirdls_data(iDay,Settings.InDir);
  if HirdlsData.Error ~= 0; 
    %problem getting data - skip
    disp('Problem getting data')
    clear HirdlsData
    continue
  end
  
  %extract prs scale for later use - it's constant by definition across the
  %whole HIRDLS mission
  Prs = HirdlsData.Prs(1,:);
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %travel angle!
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

  %find the travel direction at each point
  %assume one height level - should be pretty close
  LatScale = HirdlsData.Lat(:,20);
  LonScale = HirdlsData.Lon(:,20);
  
  %get the next point
  NextLat = circshift(LatScale,[0,-1]);
  NextLon = circshift(LonScale,[0,-1]);

  %this introduce one bad point - kill it
  LatScale = LatScale(2:end); LonScale = LonScale(2:end);
  NextLat  = NextLat( 2:end); NextLon  = NextLon( 2:end);
  
  %then find the travel direction at each point
  Azimuth = azimuth(LatScale,LonScale,NextLat,NextLon)'; 
  Azimuth = [NaN,Azimuth];
  Azimuth(1) = Azimuth(2); %should be very close

  
  %HIRDLS looks out at an angle of 47 degrees off-track behind the satellite
  %so, find the travel angle and then convert
  ViewAngleH = wrapTo180(360-((360-Azimuth)+180+47)); %degrees c/w from N
  warning off %warning constantly popping up about a future version change - supporess for now
  ViewAngleH = repmat(ViewAngleH',[1,size(Prs)]);
  warning on
  
  clear NextLat NextLon ThisLat ThisLon Azimuth BasisLevel Drift
  clear NewLon NewLat

  %z viewing angle is zero - instrument is viewing from the limb
  ViewAngleZ = zeros(size(ViewAngleH));  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %get the data we want
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
  
  
  %reduce down to just the pressure range we want
  InPrsRange = find(Prs >= min(Settings.PrsRange) ...
                  & Prs <= max(Settings.PrsRange));
  Prs = Prs(  InPrsRange);
  Lat = HirdlsData.Lat(:,InPrsRange);
  Lon = HirdlsData.Lon(:,InPrsRange);
  ViewAngleH = ViewAngleH(:,InPrsRange);
  ViewAngleZ = ViewAngleZ(:,InPrsRange);  
  clear InPrsRange
  
  %generate 30km lat array, for trimming later
  [~,idx] = min(abs(p2h(30)-Prs));
  LatScale = Lat(:,idx);
  clear idx
  
  %merge into point x height grid
  [Prs,Time] = meshgrid(Prs,HirdlsData.Time(:,1)); %we assumed time constant across a profile   earlier

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
  ShapeX  = 200; %200km along-LOS
  ShapeY  = 20;  %20km across-LOS
  ShapeZ1 = 1;   %1km vertical width in s'sphere
  ShapeZ2 = 2;   %2km vertical width in m'sphere
  
  %the values are therefore
  Weight.X = single(ones(size(Recon.x)).*ShapeX./4); %/2 for stdevs, additional /2 for width rather than distance from centre 
  Weight.Y = single(ones(size(Recon.x)).*ShapeY./4);   

  %separate z in m'sphere and s'sphere
  
  Weight.Z = single(ones(size(Recon.x))); 
  
  LowZ = find(p2h(Prs(:)) <= 60);
  Weight.Z(LowZ) = ShapeZ1./4;
  
  HighZ = find(p2h(Prs(:)) > 60);
  Weight.Z(HighZ) = ShapeZ2./4;
  
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
  clear Track DayFile
  disp([datestr(iDay),' complete'])
end
clear iDay
