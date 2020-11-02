function Model = load_jra55c(DayNumber)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load JRA55C data for a particular day, and put into common analysis format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%possible errors:
%0: success
%1. date not in valid range
%2. file not found

%get core variables - needed for model data path
CoreVars = sampling_core_variables;
CoreVars.JRA55C.Path = [LocalDataDir,'/JRA55C/'];

%settings
Settings.DataDir = CoreVars.JRA55C.Path;
Settings.DateFormat = 'mm/dd/yyyy (HH:MM)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%identify the files for this day
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


LastFile = '';
for iDay=1:1:2; %loop over the day in question and the day following, and find all the files we need
  
  %%JRA55 file format:
  %one file per ten days, all in same file except p which is level-based
  [y,m,d] = datevec(DayNumber+iDay-1);

  %first day is the day-of-month rounded down to the nearest ten, plus one
  %unless it's the 31st in a month with 31 days, then it's the 21st
  if d ~= 31;
    dd = 1+floor((d-1)/10)*10;
    StartDate = [sprintf('%04d',y),sprintf('%02d',m),sprintf('%02d',dd),'00'];
    clear dd    
    else
    StartDate = [sprintf('%04d',y),sprintf('%02d',m),'2100'];
  end

  
  %last day is the smaller of either this day plus ten or the last of the month
  %but it may shift up one if that's the end of the month
  dd = ceil(d/10)*10;
  NMonthDays = cjw_nmonthdays(m,y);
  if dd > NMonthDays; dd = NMonthDays; end
  if dd == 30 && NMonthDays == 31; dd = 31; end
  EndDate = [sprintf('%04d',y),sprintf('%02d',m),sprintf('%02d',dd),'18'];
  clear dd NMonthDays
  
  
  %hence:
  FileName.T{iDay} = ['anl_mdl.C.011_tmp.reg_tl319.',StartDate,'_',EndDate,'.nc4'];
  
end

%drop down to a single filename if the two days are in the same file
if strcmp(FileName.T{1},FileName.T{2}) == 1; NToDo = 1; else NToDo = 2; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iFile=1:1:NToDo;
  
  %retain the data form the last file
  if iFile == 2; DataSoFar = Data; end

  %identify filename
  FilePath.T = [Settings.DataDir,'/',FileName.T{iFile}];
  
  
  %get the data (or not)
  if exist(FilePath.T) ~= 0;
    %load the file
    Data.longitude = full(ncArray(FilePath.T,'g4_lon_3'));
    Data.latitude  = full(ncArray(FilePath.T,'g4_lat_2'));
    Data.level     = full(ncArray(FilePath.T,'lv_HYBL1'));
    Data.time      = full(ncArray(FilePath.T,'initial_time0'));
    Data.T         = full(ncArray(FilePath.T,'TMP_GDS4_HYBL'));
    
    %convert time to Matlab time
    t = NaN(size(Data.time,2),1);
    
    for it=1:1:numel(t); t(it) = datenum(Data.time(:,it)',Settings.DateFormat); end
    Data.time = t; clear t it
    %latitudes are monotonically decreasing, so flip in that dimension
    [~,Order] = sort(Data.latitude,'ascend');
    Data.latitude = Data.latitude(Order);
    Data.T = Data.T(:,Order,:,:);
    clear Order
    %change longitudes to run from -180 to +180
    Data.longitude(Data.longitude > 180) = Data.longitude(Data.longitude > 180)-360;
    [~,Order] = sort(Data.longitude,'ascend');
    Data.longitude = Data.longitude(Order);
    Data.T = Data.T(Order,:,:,:);
    clear Order FilePath
    
  end
  
  
  %merge the two files
  if iFile == 2;
    DataSoFar.T    = cat(4,DataSoFar.T,Data.T);
    DataSoFar.time = [DataSoFar.time;Data.time];
    Data = DataSoFar; clear DataSoFar
  end
  


end



if ~exist('Data');
  Model.Error = 2;
  return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%reformat for analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%output format:
%struct called Model
%containing fields:
%Lon   - 1d. Runs from -180 to +180
%Lat   - 1d
%Time  - 1d
%Prs   - 1d
%T     - 4d, time x lon x lat x pressure


%first, compute an average pressure scale
%we're in the stratosphere, and have much bigger errors than this elsewhere
PrsScale = squeeze(jra55_prs(10,60))';

%I'd love to keep the whole ten-day period and not do multiple load steps
%but this keeps causing out-of-memory crashes on the HPC cluster
%so drop it down to a single day
OnThisDay = find(floor(Data.time) == DayNumber);
OnThisDay = [OnThisDay;OnThisDay(end)+1]; %get the first timestep of the next day along
Data.time = Data.time(OnThisDay);
Data.T    = Data.T(:,:,:,OnThisDay);

%then, stick stuff in a struct
Model.Lon  = Data.longitude;
Model.Lat  = Data.latitude;
Model.Time = Data.time;
Model.T    = permute(Data.T,[4,1,2,3]);
Model.Prs  = PrsScale;


%success!
Model.Error = 0;
return


