function [CosmicData,OldFile] = extract_cosmic_data(MatlabDay,DataDir,OldFile,NoFlatten)

if ~exist('NoFlatten'); NoFlatten = 0; end
if ~exist('OldFile'); OldFile.Name = ''; end

CoreVars = sampling_core_variables;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%extract COSMIC data for a given day
%only runs if calling a new file to previous iteration
%
%Corwin Wright, corwin.wright@trinity.oxon.org
%27/DEC/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%errors:
%0. no error - successful
%1. no file found
%2. no data on this day
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find the file (if it's not the one from last time) and load it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfield(OldFile,'Name'); OldFile.Name = ' '; end; %always load file afresh if old file not specified

%identify file for this day
[y,~,~]   = datevec(MatlabDay);
dayno = date2doy(MatlabDay);

FileString = strcat('cosmic_',sprintf('%04d',y),'_',sprintf('%03d',dayno),'.mat');
if strcmp(FileString,OldFile.Name) == 0;
  %new file - load it up
  
  FileName = wildcardsearch([DataDir,'/',sprintf('%04d',y)],FileString);
  
  if numel(FileName) == 0;  %no file
    CosmicData.Error = 1;
    return;
  end;
  
  %file exists! load data
  AllCosmicData = load(FileName{1}); 
  AllCosmicData = AllCosmicData.Data;
  clear FileName;
  
  OldFile.Name = FileString;
  
else
  %same as last call - don't reload
  AllCosmicData = OldFile.Data;
end
clear m y DataDir FileString

%keep for next time around
OldFile.Data = AllCosmicData;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%interpolate onto common height scale
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Data = AllCosmicData;

%COSMIC files have a massive number of levels, but in the stratosphere
%resolution is 1.5km due to optics. so, to save sampling time, 
%interpolate to constant height scale
%we'll use a 1/16 decade scale (about 1km) - this is still oversampled relative
%to resolution
Settings.LogPSpacing = 1/16;

%COSMIC data can be pretty NaNny for methodological reasons
%maximum frac of NaNs in profile:
Settings.MaxNaNFrac = 0.75;%



%interpolate all profiles onto a constant log-pressure scale
Settings.PrsRange   = CoreVars.Cosmic.HeightRange;
PrsScale = 10.^(log10(Settings.PrsRange(2)):Settings.LogPSpacing:log10(Settings.PrsRange(1)));
Profs.Prs = repmat(PrsScale',1,size(Data,2));


VarsIn = {'Lon','Lat','Time','Temp'};
VarsOut = {'Lon','Lat','Time','T'};
for iVar=1:1:numel(VarsIn);
  Grid = Profs.Prs.*NaN;
  for iProfile=1:1:size(Data,2);
    
    if Data(iProfile).QC ~= 0; continue; end %failed quality control at CDAAC
    
    %make sure the scales are monotonically increasing by...
    %ordering
    x = Data(iProfile).Pres;
    v = Data(iProfile).(VarsIn{iVar});
    [~,order] = sort(x,'ascend');
    x = x(order); v = v(order);
    %and de-duplicating (this method removes one point at the end, but we're
    %interpolating ~3000 onto ~45 so that's not a big issue)
    NotDupe = find(diff(x) ~= 0);
    x = x(NotDupe); v = v(NotDupe);
    
    %remove NaNs, and skip if above a cutoff
    Good = find(~isnan(x+v));
    Frac = numel(Good)./numel(x);
    if Frac < Settings.MaxNaNFrac; continue; end
    x = x(Good); v = v(Good);
    
    Grid(:,iProfile) = interp1(x,v,PrsScale);
  end
  Profs.(VarsOut{iVar}) = Grid;
  clear Grid iVar Frac x v order NotDupe Good
end; clear iProfile

CosmicData = Profs;


if NoFlatten ==0;
  %produce lat and lon scale
  CosmicData.Lat = nanmean(CosmicData.Lat,1);
  CosmicData.Lon = nanmean(CosmicData.Lon,1);
end


%flip it all, to match other routines
CosmicData.Lon  = CosmicData.Lon';
CosmicData.Lat  = CosmicData.Lat';
CosmicData.Time = CosmicData.Time';
CosmicData.Prs  = CosmicData.Prs';
CosmicData.T    = CosmicData.T'+273.15; %convert from C to K at same time




CosmicData.Error = 0; %it worked!




return
end
