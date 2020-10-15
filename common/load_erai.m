function Model = load_erai(DayNumber)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load ERAI data for a particular day, and put into common analysis format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%possible errors:
%0: success
%1. date not in valid range
%2. file not found

%get core variables - needed for model data path
CoreVars = sampling_core_variables;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load the data for this day
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%ERAi file format:
%monthly files, with p based on a standard grid
Settings.DataDir    = CoreVars.ERAI.Path;

%load the appropriate files
[y,m,~] = datevec(DayNumber);

FileName.T = ['era-interim_an_ml_',sprintf('%04d',y),'-',sprintf('%02d',m),'_T.nc'];
FileName.P = ['era-interim_an_ml1_',sprintf('%04d',y),'-',sprintf('%02d',m),'_lnsp.nc'];

  
%try and load the file
FilePath.T = [Settings.DataDir,'/',FileName.T];
FilePath.P = [Settings.DataDir,'/',FileName.P];


if exist(FilePath.T) ~= 0 & exist(FilePath.P) ~= 0;
  %load the file
  Data.longitude = full(ncArray(FilePath.T,'longitude'));
  Data.latitude  = full(ncArray(FilePath.T,'latitude'));
  Data.level     = full(ncArray(FilePath.T,'level'));
  Data.time      = full(ncArray(FilePath.T,'time'));
  Data.T         = full(ncArray(FilePath.T,'t'));
  Data.lnsp      = full(ncArray(FilePath.P,'lnsp'));
  
  %convert time to Matlab time
  %default format is hours since 1900-01-01 00:00:00
  Data.time = datenum(1900,1,1,double(Data.time),0,0);
  
  %latitudes are monotonically decreasing, so flip in that dimension
  [~,Order] = sort(Data.latitude,'ascend');
  Data.latitude = Data.latitude(Order);
  Data.T = Data.T(:,Order,:,:);
  Data.lnsp = Data.lnsp(:,Order,:);
  clear Order
  
  %change longitudes to run from -180 to +180
  Data.longitude(Data.longitude > 180) = Data.longitude(Data.longitude > 180)-360;
  [~,Order] = sort(Data.longitude,'ascend');
  Data.longitude = Data.longitude(Order);
  Data.T = Data.T(Order,:,:,:);
  Data.lnsp = Data.lnsp(Order,:,:);
  clear Order
  
  clear FileName y  FilePath
end


%if we need to load another month to get the next day, load it, process it, and
%cat it to the data
[y2,m2,~] = datevec(DayNumber+1);
if m2 ~= m;

  
  FileName.T = ['era-interim_an_ml_',sprintf('%04d',y2),'-',sprintf('%02d',m2),'_T.nc'];
  FileName.P = ['era-interim_an_ml1_',sprintf('%04d',y2),'-',sprintf('%02d',m2),'_lnsp.nc'];
  FilePath.T = [Settings.DataDir,'/',FileName.T];
  FilePath.P = [Settings.DataDir,'/',FileName.P];
  
  if exist(FilePath.T) ~= 0 & exist(FilePath.P) ~= 0;  %if the file doesn't exist, just stick with what we have - it's most of a day anyway!

    DataSoFar = Data; clear Data;
    %load the file
    Data.longitude = full(ncArray(FilePath.T,'longitude'));
    Data.latitude  = full(ncArray(FilePath.T,'latitude'));
    Data.level     = full(ncArray(FilePath.T,'level'));
    Data.time      = full(ncArray(FilePath.T,'time'));
    Data.T         = full(ncArray(FilePath.T,'t'));
    Data.lnsp      = full(ncArray(FilePath.P,'lnsp'));
    
    %convert time to Matlab time
    %default format is hours since 1900-01-01 00:00:00
    Data.time = datenum(1900,1,1,double(Data.time),0,0);
    
    %latitudes are monotonically decreasing, so flip in that dimension
    [~,Order] = sort(Data.latitude,'ascend');
    Data.latitude = Data.latitude(Order);
    Data.T = Data.T(:,Order,:,:);
    Data.lnsp = Data.lnsp(:,Order,:);
    clear Order
    
    %change longitudes to run from -180 to +180
    Data.longitude(Data.longitude > 180) = Data.longitude(Data.longitude > 180)-360;
    [~,Order] = sort(Data.longitude,'ascend');
    Data.longitude = Data.longitude(Order);
    Data.T = Data.T(Order,:,:,:);
    Data.lnsp = Data.lnsp(Order,:,:);
    clear Order
    
    clear FileName y m  FilePath
    
    %cat it
    DataSoFar.time = [DataSoFar.time;Data.time];
    DataSoFar.T    = cat(4,DataSoFar.T,Data.T);
    
    %done
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
%we're in the stratosphere, and have much bigger errors than thuis elsewhere
PrsScale = ecmwf_prs_v2(squeeze(Data.lnsp(1,1,1)),60)';

%the files we are using are not complete, to save filespace. trim down to what we actually have
PrsScale = PrsScale(1:1:32);
PrsScale(1) = 0.1; %because it is, but routine returns a NaN

%I'd love to keep the whole month and not do multiple load steps
%but this keeps causing out-of-memory crashes on the HPC cluster
%so drop it down to a single day, including the first timestpe of the next day
OnThisDay = find(floor(Data.time) == DayNumber);
OnThisDay = [OnThisDay',max(OnThisDay)+1]; %first timestep of next day
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


