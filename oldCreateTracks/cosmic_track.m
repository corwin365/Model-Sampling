clearvars
addpath(genpath('../'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%produce OIF data for COSMIC
%
%Corwin Wright, c.wright@bath.ac.uk, 2023/02/03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%dataset identifier
Settings.Instrument = 'COSMIC';

%where do the input files live?
Settings.InDir = [LocalDataDir,'/COSMIC/daily_atmPrf'];

%geolocation - which data should we include?
%for all except HeightRange, we include any wholegranule including these
%for HeightRange, we will trim the granules in height to just this range
Settings.LatRange    = [-90,90];
Settings.LonRange    = [-180,180];
Settings.TimeRange   = datenum(2020,1,24);
Settings.HeightRange = [0,50]; %km

%path handling internal to routine
[~,CoreSettings] = sampling_core_v2(' ',' ',0,'GetSettings',true);
Settings.OutDir  = [CoreSettings.MasterPath,'/tracks/',Settings.Instrument,'/'];
clear CoreSettings

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%extract lat,lon,z of each COSMIC point, to allow model sampling
%store in daily files of a common format
%
%unlike the others, also store T, since it's downsampled.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%COSMIC files have a massive number of levels, but in the stratosphere
%resolution is 1.5km due to optics. so, to save sampling time, 
%interpolate to constant height scale
%we'll use a 1/32 decade scale (about 0.5km) - this is still oversampled relative
%to resolution
Settings.LogPSpacing = 1/32;

%COSMIC data can be pretty NaNny for methodological reasons
%maximum frac of NaNs in profile:
Settings.MaxNaNFrac = 0.75;%

for iDay=min(Settings.TimeRange):1:max(Settings.TimeRange);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %generate file name
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  %create a directory to put the data in, if it doesn't exist
  if exist(Settings.OutDir,'dir') ~= 7; mkdir(Settings.OutDir); end

  %and name the output file
  DayFile = [Settings.OutDir,'track_',Settings.Instrument,'_',num2str(iDay),'.mat'];
  if exist(DayFile); continue; end
  
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
  PrsScale = 10.^(log10(h2p(max(Settings.HeightRange))):Settings.LogPSpacing:log10(h2p(min(Settings.HeightRange))));
  Profs.Prs = repmat(PrsScale',1,size(Data,2));

  VarsIn = {'Lon','Lat','Azim','Time','Temp'};
  VarsOut = {'Lon','Lat','Az','Time','T'};
  for iVar=numel(VarsIn):-1:1;
    Grid = NaN(size(Profs.Prs));
    for iProfile=1:1:size(Data,2);

      %if time, duplicate out to a full profile as we only stored one point in these files
      if strcmp(VarsIn{iVar},'Time');
        t = Data(iProfile).(VarsIn{iVar});
        Data(iProfile).(VarsIn{iVar}) = ones(size(Data(iProfile).Lat)).*t;
      end
      
      %make sure the scales are monotonically increasing by...
      %ordering
      x = Data(iProfile).Pres;
      v = Data(iProfile).(VarsIn{iVar});
      [~,order] = sort(x,'ascend');
      x = x(order); v = v(order);
      %and de-duplicating (this method removes one point at the end, but we're
      %interpolating ~3000 onto ~45 so that's not a big issue)
      NotDupe = find(abs(diff(x)) > 0);
      x = x(NotDupe); v = v(NotDupe);
      
      %remove NaNs, and skip profile if above a cutoff fraction
      Good = find(~isnan(x+v) > 0);
      Frac = numel(Good)./numel(x);
      if Frac < Settings.MaxNaNFrac | numel(Good) < 2;     continue; end
      x = x(Good); v = v(Good);


      Grid(:,iProfile) = interp1(x,v,PrsScale);
     
    end

    Profs.(VarsOut{iVar}) = Grid;
    clear Grid iVar Frac x v order NotDupe Good
  end; clear iProfile

  %azimuths are c/w from N. which is what we want. yay!

% % % %   %store QC flag, for later use
% % % %   %keep both good and bad data as we may want to compare perfect-QC COSMIC with
% % % %   %failed-QC COSMIC at a later date. remove them at the analysis stage if not.
% % % %   Profs.QC = NaN(size(Profs.Prs));
% % % %   for iQC=1:1:size(Profs.Prs,2);
% % % %     Profs.QC(:,iQC) = Data(iQC).QC;
% % % %   end
% % % %   VarsOut = cat(2,VarsOut,'QC');
  
  
  %tidy up
  clear dn PrsScale VarsIn y InFile

  %remove NaN profiles
  NaNProf = sum(Profs.T,1) + sum(Profs.Lon,1) + sum(Profs.Lat,1) + sum(Profs.Az,1) + sum(Profs.Time,1); %can omit prs as we made it
  Good = find(~isnan(NaNProf));  
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
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
   
  
  %put data into a line
  Prs  = Prs(:);
  Time = Time(:);
  Lat  = Lat(:);
  Lon  = Lon(:);  
  T    = T(:);
% % %   QC   = QC(:);
  ViewAngleH = ViewAngleH(:);
  ViewAngleZ = ViewAngleZ(:);  
  
  %then into an array
  Track.T          = single(T);
  Track.Lat        = single(Lat);
  Track.Lon        = single(Lon);
  Track.Prs        = single(Prs);
% % %   Track.QC         = single(QC);
  Track.Time       = Time; %needs to be double
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
