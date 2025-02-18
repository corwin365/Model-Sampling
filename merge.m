 clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%merge ERA5-as-AIRS files sampled for Neil into daily files
%
%Corwin wright, c.wright@bath.ac.uk
%2024/01/28
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Settings.DayScale = datenum(2002,12,8):1:datenum(2002,12,14);
Settings.DataDir  = [LocalDataDir,'/corwin/sampling_project/output/AIRS3D/era5_levante/'];
Settings.OutDir   = [LocalDataDir,'/corwin/sampling_project/output/AIRS3D/era5_merged/'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% loop over days, merging found files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iDay=1:1:numel(Settings.DayScale)
%  try  

  % %get list of files
  % Files = wildcardsearch(Settings.DataDir,['*',num2str(Settings.DayScale(iDay)),'*']);
  % if numel(Files) == 0; disp([datestr(Settings.DayScale(iDay)),' skipped, no data']); clear Files; continue; end
  
  %name outfile
  OutFile = [Settings.OutDir,'/',yyyyDDD(Settings.DayScale(iDay)),'.mat'];
  if exist(OutFile,'file'); disp(['Already done ',datestr(Settings.DayScale(iDay))]); continue; end

  %create results array for this day

%    T       = NaN(240,135,90,21,'single');
%    Z       = [12,15,18,21,24,27,30,33,36,39,42,45,48,51,54,57,60,65,70,75,80];
  
  T       = NaN(240,135,90,14,'single');
  Z       = [21,23,27,30,33,36,39,42,45,48,51,54,57,60];
  TSimple = T;
  
  %loop over files, storing in array
  % for iFile=1:1:numel(Files)
  for iG=1:1:240;
  
    % %get granule number from file name
    % G = Files{iFile}; G = str2num(G(end-6:end-4));
    %generate filename
    File = [Settings.DataDir,'sampled_',num2str(Settings.DayScale(iDay)),'_subset',sprintf('%06d',iG),'.mat'];

    %load file and see if the format is full with metadata or minimal T-only
    try; Data = load(File); catch; clear File; continue; end
    if isfield(Data,'Sampled_Data')
      T(      iG,:,:,:) = permute(Data.Sampled_Data.T,[2,3,1]);
      TSimple(iG,:,:,:) = permute(Data.Sampled_Data.TSimple,[2,3,1]);
    else
      T(      iG,:,:,:) = permute(Data.T,      [2,3,1]);
      TSimple(iG,:,:,:) = permute(Data.TSimple,[2,3,1]);
    end
    


    %delete the original file
    delete(File)

    %tidy up and loop
    clear Data File    

  end; clear iG

  %save
  save(OutFile,'Z','T','TSimple');
  disp([datestr(Settings.DayScale(iDay)),' saved'])


%  catch;  disp(['Error with levels on ',datestr(Settings.DayScale(iDay))]);
%  end

end; clear iDay
