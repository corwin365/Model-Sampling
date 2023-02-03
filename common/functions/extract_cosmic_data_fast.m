function [Cosmic,OldFile] = extract_cosmic_data(MatlabDay,NoFlatten)

if ~exist('NoFlatten'); NoFlatten = 0; end
if ~exist('OldFile'); OldFile.Name = ''; end

CoreVars = sampling_core_variables;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%extract COSMIC data for a given day
%
%uses pre-interpolated data in the track files to speed up loading, as raw
%COSMIC data is very large and slow to load
%
%Corwin Wright, corwin.wright@trinity.oxon.org
%16/JAN/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%errors:
%0. no error - successful
%1. no file found



%find the file (if it's not the one from last time) and load it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TrackDir = [CoreVars.MasterPath,'/tracks/COSMIC'];
DayFile = [TrackDir,'/track_cosmic_',num2str(MatlabDay),'.mat'];
if ~exist(DayFile); Cosmic.Error = 1;  return; end


%load and format data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load
%%%%%%%%%%%%%%%%%%%%%%%
Data = load(DayFile);

%reshape
%%%%%%%%%%%%%%%%%%%%%%%

%first, work out the maximum profiles per day
NProfs  = max(Data.Recon.x);
%and the maximum number of levels
NLevels = max(Data.Recon.z);

%hence, create arrays to store the reconstructed data
Cosmic         = struct();
Cosmic.T       = NaN(NProfs,NLevels);
Cosmic.Lat     = Cosmic.T;
Cosmic.Lon     = Cosmic.T;
Cosmic.Time    = Cosmic.T;
Cosmic.Prs     = Cosmic.T;
Cosmic.QC      = Cosmic.T;

%% extract


for iPoint=1:1:numel(Data.Track.Lat);
  
  x = Data.Recon.x(iPoint);
  z = Data.Recon.z(iPoint);
  
  Cosmic.T(      x,z) = Data.Track.T(      iPoint)+273.15;
  Cosmic.Lat(    x,z) = Data.Track.Lat(    iPoint);
  Cosmic.Lon(    x,z) = Data.Track.Lon(    iPoint);
  Cosmic.Time(   x,z) = Data.Track.Time(   iPoint);
  Cosmic.Prs(    x,z) = Data.Track.Prs(    iPoint); %we worked in logs
  Cosmic.QC(     x,z) = Data.Track.QC(     iPoint);
  
end

if NoFlatten ==0;
  %produce lat and lon scale
  Cosmic.Lat = nanmean(Cosmic.Lat,1);
  Cosmic.Lon = nanmean(Cosmic.Lon,1);
  Cosmic.QC  = nanmean(Cosmic.QC, 1);
end


%apply QC flags
Cosmic.T(Cosmic.QC == 1) = NaN;

Cosmic.Error = 0;

return
end
