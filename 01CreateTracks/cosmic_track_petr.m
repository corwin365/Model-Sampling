clearvars
addpath('../common');
CoreVars = sampling_core_variables;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%extract lat,lon,z of each COSMIC point, to allow model sampling
%store in daily files of a common format
%
%This is a modified version for just Andes/Peninsula, but at full res.
%
%
%unlike the others, also store T, since it's downsampled.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%general settings
Settings.Instrument = 'cosmic_petr';
Settings.InDir      = CoreVars.Cosmic.Path;
Settings.OutDir     = [CoreVars.MasterPath,'/tracks/COSMICpetr/'];
Settings.PrsRange   = CoreVars.Cosmic.HeightRange;
Settings.TimeRange  = datenum(2010,10,8):1:datenum(2010,10,19)

%restrict to study region
Settings.LatRange   = [-75,-40];
Settings.LonRange   = [-85,-35];

%interpolate to constant height scale
%we'll use a 1/1000 decade scale - this is higher than any actual data in the file (just)
Settings.LogPSpacing = 1/1000;

%COSMIC data can be pretty NaNny for methodological reasons
%maximum frac of NaNs in profile:
Settings.MaxNaNFrac = 0.75;%

for jDay=1:1:numel(Settings.TimeRange);
  
  iDay = Settings.TimeRange(jDay);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% import data
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  %identify file
  [y,~,~] = datevec(iDay);
  dn = date2doy(iDay);
  InFile = [Settings.InDir,'/',sprintf('%04d',y),'/cosmic_',sprintf('%04d',y),'_',sprintf('%03d',dn),'.mat'];

  if ~exist(InFile); continue; end
  load(InFile)



  %interpolate all profiles onto a constant log-pressure scale
  PrsScale = 10.^(log10(Settings.PrsRange(2)):Settings.LogPSpacing:log10(Settings.PrsRange(1)));
  Profs.Prs = repmat(PrsScale',1,size(Data,2)); 

  VarsIn  = {'Lon','Lat','Azim','Time','Temp'};
  VarsOut = {'Lon','Lat','Az','Time','T'};
  for iVar=numel(VarsIn):-1:1;
    Grid = NaN(size(Profs.Prs));
    for iProfile=1:1:size(Data,2);
      
      %make sure the scales are monotonically increasing by...
      %ordering
      x = Data(iProfile).Pres;
      v = Data(iProfile).(VarsIn{iVar});
      [~,order] = sort(x,'ascend');
      x = x(order); v = v(order);
      
      %and de-duplicating (this method removes one point at the end, not a big issue)
      NotDupe = find(abs(diff(x)) > 0);
      x = x(NotDupe); v = v(NotDupe);
      
      %remove NaNs, and skip profile if above a cutoff fraction
      Good = find(~isnan(x+v) > 0);
      Frac = numel(Good)./numel(x);
      if Frac < Settings.MaxNaNFrac; continue; end
      x = x(Good); v = v(Good);

      Grid(:,iProfile) = interp1(double(x),v,PrsScale);
     
    end

    Profs.(VarsOut{iVar}) = Grid;
    clear Grid iVar Frac x v order NotDupe Good
  end; clear iProfile

  %azimuths are c/w from N. which is what we want. yay!

  %store QC flag, for later use
  %keep both good and bad data as we may want to compare perfect-QC COSMIC with
  %failed-QC COSMIC at a later date. remove them at the analysis stage if not.
  Profs.QC = NaN(size(Profs.Prs));
  for iQC=1:1:size(Profs.Prs,2);
    Profs.QC(:,iQC) = Data(iQC).QC;
  end
  VarsOut = cat(2,VarsOut,'QC');
  
  
  %tidy up
  clear dn PrsScale VarsIn y InFile

  %remove NaN profiles, and discard any profiles not entirely in the study region
  NaNProf = sum(Profs.T,1) + sum(Profs.Lon,1) + sum(Profs.Lat,1) + sum(Profs.Az,1) + sum(Profs.Time,1); %can omit prs as we made it
  Good = find(~isnan(NaNProf));  
  
  InRegion = [];
  for iProf=1:1:size(Profs.Lat,2)
    if nanmin(Profs.Lon(:,iProf)) > min(Settings.LonRange) ...
     & nanmin(Profs.Lat(:,iProf)) > min(Settings.LatRange) ...
     & nanmax(Profs.Lon(:,iProf)) < max(Settings.LonRange) ...
     & nanmax(Profs.Lat(:,iProf)) < max(Settings.LatRange);
      InRegion(end+1) = iProf;
    end
  end
  
  
  Good = intersect(Good,InRegion);
  Profs.Prs = Profs.Prs(:,Good); %handled separately
  for iVar=1:1:numel(VarsOut);
    Var = Profs.(VarsOut{iVar});
    Var = Var(:,Good);
    Profs.(VarsOut{iVar}) = Var;
  end; clear iVar Var NaNProf Good VarsOut
  
 
  
  

  %and rename/reshape variables to reuse old code below (copied from HIRDLS routine)
  Lat  = Profs.Lat';
  Lon  = Profs.Lon';
  Time = Profs.Time';
  T    = Profs.T';
  ViewAngleH = Profs.Az';
  Prs  = Profs.Prs';
  QC   = Profs.QC';
  
  %viewing angle from above isn't terribly meaningful for COSMIC
  %so set it to zero
  ViewAngleZ = ViewAngleH; ViewAngleZ(:) = 0;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %prepare the necessary information to reconstruct the granules from 1D data
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

  x = 1:1:size(Lat,1);  
  z = 1:1:size(Lat,2);
  [z,x] = meshgrid(z,x);
  Recon.x = single(x(:));
  Recon.z = single(z(:));
  clear x z

  
  %also store a unique average time for each profile, so that when we 
  %split the data into hourly chunks none span these time interfaces
  t = repmat((nanmean(Profs.Time,1))',1,size(Time,2));
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %define the approximate weighting blob
  %numbers I'm using are ~2 stdev, so divide appropriately
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
  
  %approximate shape of blob
  ShapeX = 270; %270km along-LOS
  ShapeY = 1.5; %1.5km across-LOS
  ShapeZ = 1.5; %1.5km vertical width
  
  %the values are therefore
  Weight.X = single(ones(size(Recon.x)).*ShapeX./4); %/2 for stdevs, additional /2 for width rather than distance from centre 
  Weight.Y = single(ones(size(Recon.x)).*ShapeY./4);   
  Weight.Z = single(ones(size(Recon.x)).*ShapeZ./4);     
  

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %save!
  %
  %split into one-hour chunks to keep memory manageable
  %for the big high-res runs Chris and Annelize provided
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
   
  
  %put data into a line
  Prs  = Prs(:);
  Time = Time(:);
  Lat  = Lat(:);
  Lon  = Lon(:);  
  T    = T(:);
  QC   = QC(:);
  ViewAngleH = ViewAngleH(:);
  ViewAngleZ = ViewAngleZ(:);  
  
  %then into an array
  Track.T          = single(T);
  Track.Lat        = single(Lat);
  Track.Lon        = single(Lon);
  Track.Prs        = single(Prs);
  Track.QC         = single(QC);
  Track.Time       = Time; %needs to be double
  Track.ViewAngleH = single(ViewAngleH);
  Track.ViewAngleZ = single(ViewAngleZ);
  
  clear Lat Lon Prs Time

  %and save it
  TrackStore  = Track;  clear Track
  ReconStore  = Recon;  clear Recon
  WeightStore = Weight; clear Weight
  t = t(:);
  [~,~,~,h,~,~] = datevec(t);
  for iChunk=1:1:24;
    
    %select all the points in this hour
    InHour = find(h == iChunk);
    
    %create track and recon arrays considering of these points
    Track  = reduce_struct( TrackStore,InHour);
    Recon  = reduce_struct( ReconStore,InHour);
    Weight = reduce_struct(WeightStore,InHour);
    
    %and store
    OutFile = [Settings.OutDir,'track_',Settings.Instrument,'_',num2str(iDay),'_g',sprintf('%03d',iChunk),'.mat'];
    save(OutFile,'Track','Recon','Weight');
  end
  
  
  
  
  
  %tidy up, then done
  clear Track DayFile
  disp([datestr(iDay),' complete'])

end
clear iDay
