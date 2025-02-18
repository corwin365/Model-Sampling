clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%generate lists of:
%   - AIRS trackfiles
%   - ERA5 T files
%  needed to fill gaps in the climatology generated so far
%
%Corwin Wright, c.wright@bath.ac.uk
%2024/MAR/02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%path to data store
Settings.OutputDir = [LocalDataDir,'/corwin/sampling_project/output/AIRS3D/era5_merged/test'];

%path to sampler output
Settings.SampleDir  = [LocalDataDir,'/corwin/sampling_project/output/AIRS3D/era5_levante/'];

%list of granules to operate on
Settings.GranuleList =  'airs_missing.mat';

%prepare access to the sampler
addpath(genpath('../'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get list of unique days
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

List = load(Settings.GranuleList); List = List.Needed;

Days = unique(List(:,1));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% process!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iDay=1:1:Days

  %% load the original saved day we're plugging the gaps in
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %load the storage file for this day
  StoreFile = [Settings.OutputDir ,'/',yyyyDDD(Days(iDay)),'.mat'];
  if ~exist(StoreFile,'file')
    %all granules are 'missing' - create a new file
    NewFile = 1;
  else
    %load the file and identify the missing granules
    try
      StoreFile = load(StoreFile);
      NewFile = 0;
    catch; 
      NewFile = 1;
    end
  end
  
  if NewFile == 1; continue; end %address this later to fill the fully-empty days


  %% find missing granules on this day
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  Granules = List(Days(iDay) == List(:,1) ,2);

 
  %% sample them
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  for iG=1:1:numel(Granules)
    cd ..
    sampling_core_v3('AIRS3D','era5_levante',Days(iDay),'SubSet',Granules(iG),'Parallel',true,'SaveTOnly',true);
    cd gapfill
  end



  if NewFile == 0;
   
  
    %% if we have the day already, plug the gaps
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for iG=1:1:numel(Granules)
      %load the granule
      GG = load([Settings.SampleDir,'sampled_',num2str(Days(iDay)),'_subset',sprintf('%06d',Granules(iG)),'.mat']);

      %stick the new data in
      StoreFile.T(      Granules(iG),:,:,:) = permute(GG.T,      [2,3,1]);
      StoreFile.TSimple(Granules(iG),:,:,:) = permute(GG.TSimple,[2,3,1]);
      T       = StoreFile.T;
      TSimple = StoreFile.TSimple;
    end
  else
  
    %% cretae a new day and fill it
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
stop
  end
  
  
  %% save the file
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  save([Settings.OutputDir ,'/',yyyyDDD(Days(iDay)),'.mat'],'T','TSimple');
  disp([datestr(Days(iDay)),' updated and saved'])


end; clear iDay
