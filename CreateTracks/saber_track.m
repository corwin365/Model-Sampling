clearvars
addpath('../common');
CoreVars = sampling_core_variables;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%extract lat,lon,z of each SABER point, to allow model sampling
%store in daily files of a common format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%information needed. 
Settings.Instrument = 'saber';
Settings.InDir      = CoreVars.Saber.Path;
Settings.OutDir     = [CoreVars.MasterPath,'/tracks/SABER/'];
Settings.PrsRange   = CoreVars.Saber.HeightRange; 
Settings.TimeRange  = CoreVars.Saber.TimeRange; 


OldFile = '';
for iDay=Settings.TimeRange(1):1:Settings.TimeRange(2);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %generate file name
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  DayFile = [Settings.OutDir,'track_',Settings.Instrument,'_',num2str(iDay),'.mat'];
  if exist(DayFile) == 2;
    continue;
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %import data
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  [SaberData,OldFile] = extract_saber_data(iDay,Settings.InDir,OldFile,1);

  if SaberData.Error ~= 0; 
    disp(['No data found for ',datestr(iDay)])
    %problem getting data - skip
    clear SaberData
    continue
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %get the data we want
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
  
  %drop unwated vars
  SaberData = rmfield(SaberData,{'Height','Dens','T'});
  
  %tranpose the rest, to reuse old code below
  Fields = fieldnames(SaberData);
  for iField=1:1:numel(Fields)
    SaberData.(Fields{iField}) = SaberData.(Fields{iField})';
  end
  
  

  %reduce down to just the pressure range we want
  InPrsRange = find(nanmean(SaberData.Prs,2) >= min(Settings.PrsRange) ...
                  & nanmean(SaberData.Prs,2) <= max(Settings.PrsRange));
  SaberData.Prs  = SaberData.Prs( InPrsRange,:);
  SaberData.Lat  = SaberData.Lat( InPrsRange,:);
  SaberData.Lon  = SaberData.Lon( InPrsRange,:);
  SaberData.Time = SaberData.Time(InPrsRange,:);  
  clear InPrsRange

  
  %generate 30km geolocation arrays
  %we'll use these to work out the travel vector of the satellite
  [~,idx] = min(abs(p2h(30)-nanmean(SaberData.Prs,2)));
  LatScale = SaberData.Lat(idx,:);
  LonScale = SaberData.Lon(idx,:);
  clear idx
  
% %   %put prs into point x height grid
% %   [~,SaberData.Prs] = meshgrid(1:1:size(SaberData.Lat,2),SaberData.Prs);   
 
  %and rename arrays
  Lat  = SaberData.Lat';              
  Lon  = SaberData.Lon';
  Time = SaberData.Time';
  Prs  = SaberData.Prs';
  clear SaberData
  
  %horizontal viewing angle is 90 deg off track, to the right
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %get the next point
  NextLat = circshift(LatScale,[0,-1]);
  NextLon = circshift(LonScale,[0,-1]);

  %this introduce one bad point - kill it
  LatScale = LatScale(2:end); LonScale = LonScale(2:end);
  NextLat  = NextLat( 2:end); NextLon  = NextLon( 2:end);
  
  %then find the travel direction at each point
  Azimuth = azimuth(LatScale,LonScale,NextLat,NextLon); 
  Azimuth = [NaN,Azimuth];
  Azimuth(1) = Azimuth(2); %should be very close

  ViewAngleH = wrapTo180(360-((360-Azimuth)-90)); %degrees c/w from N
  ViewAngleH = repmat(ViewAngleH',1,size(Prs,2));
  
  %z viewing angle is zero - instrument is viewing from the limb
  ViewAngleZ = zeros(size(ViewAngleH));


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
  ShapeX = 300; %300km along-LOS
  ShapeY = 50;  %50km across-LOS
  ShapeZ = 2;   %2km vertical width
  
  %the values are therefore
  Weight.X = single(ones(size(Recon.x)).*ShapeX./4); %/2 for stdevs, additional /2 for width rather than distance from centre 
  Weight.Y = single(ones(size(Recon.x)).*ShapeY./4);   
  Weight.Z = single(ones(size(Recon.x)).*ShapeZ./4);   

  

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
  Track.Lat  = single(Lat);
  Track.Lon  = single(Lon);
  Track.Prs  = single(Prs);
  Track.Time = Time; %needs to be double

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