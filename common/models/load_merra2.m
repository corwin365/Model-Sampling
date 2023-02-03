function Model = load_merra2(DayNumber)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load MERRA2 data for a particular day, and put into common analysis format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%possible errors:
%0: success
%1. date not in valid range
%2. file not found

%model data path
ModelPath = [LocalDataDir,'/MERRA2/'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load the data for this day
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%day in question
[y1,m1,d1] = datevec(DayNumber);
if     y1 <= 2000; Version1 = 200;
elseif y1 <= 2010; Version1 = 300;
else               Version1 = 400;
end

%next day, to get the zero-hour timestep
[y2,m2,d2] = datevec(DayNumber+1);
if     y2 <= 2000; Version2 = 200;
elseif y2 <= 2010; Version2 = 300;
else               Version2 = 400;
end


FilePath1 = [ModelPath,'/MERRA2_',num2str(Version1),'.inst3_3d_asm_Nv.',sprintf('%04d',y1),sprintf('%02d',m1),sprintf('%02d',d1),'.SUB.nc4'];
FilePath2 = [ModelPath,'/MERRA2_',num2str(Version2),'.inst3_3d_asm_Nv.',sprintf('%04d',y2),sprintf('%02d',m2),sprintf('%02d',d2),'.SUB.nc4'];


if ~exist(FilePath1);
  Model.Error = 2; 
  return 
end
if ~exist(FilePath2);
  Model.Error = 2; 
  return 
end
 

Data1 = cjw_readnetCDF(FilePath1);
Data2 = cjw_readnetCDF(FilePath2);

Data1.time = DayNumber +     Data1.time/60/24;
Data2.time = DayNumber + 1 + Data2.time/60/24;

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

%combine the first time step of day 2
Data1.time(end+1) = Data2.time(1);
Data1.T(:,:,:,end+1) = Data2.T(:,:,:,1);
Data = Data1; clear Data1 Data2


%these are the TOPS of the levels
Model.Prs = [0.0100,0.0200,0.0327,0.0476,0.0660,0.0893,0.1197,0.1595,0.2113,0.2785,0.3650,0.4758,0.6168,0.7951,1.0194,1.3005,1.6508,2.0850,2.6202,3.2764,4.0766,5.0468,6.2168,7.6198,9.2929,11.2769,13.6434,16.4571,19.7916,23.7304,28.3678,33.8100,40.1754,47.6439,56.3879,66.6034,78.5123,92.3657,108.6630,127.8370,150.3930,176.9300,208.1520,244.8750,288.0830,337.5000,375.0000,412.5000,450.0000,487.5000,525.0000,562.5000,600.0000,637.5000,675.0000,700.0000,725.0000,750.0000,775.0000,800.0000,820.0000,835.0000,850.0000,865.0000,880.0000,895.0000,910.0000,925.0000,940.0000,955.0000,970.0000,985.0000]; 

%and so these are the level-centres
Delta = [diff(Model.Prs),15]./2;

%and so...
Model.Prs = Model.Prs+Delta;


%then, stick stuff in a struct
Model.Lon  = Data.lon;
Model.Lat  = Data.lat;
Model.Time = Data.time;
Model.T    = permute(Data.T,[4,1,2,3]);

%CFSR runs from -180 to +180, so we don't need to fix the lons

%success!
Model.Error = 0;
% return


