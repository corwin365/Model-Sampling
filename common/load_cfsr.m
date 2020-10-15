function Model = load_cfsr(DayNumber)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load CFSR data for a particular day, and put into common analysis format
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

%%CFSR file format:
%6-hourly files, with T and p in seperate files
[y1,m1,d1] = datevec(DayNumber);
FilePath1 = [CoreVars.CFSR.Path,sprintf('%04d',y1),'/',sprintf('%02d',m1),'/'];
%also need the first timestep of the next day
[y2,m2,d2] = datevec(DayNumber+1);
FilePath2 = [CoreVars.CFSR.Path,sprintf('%04d',y2),'/',sprintf('%02d',m2),'/'];

for iFile=1:1:5;
  try
    
    %load data for this timestep   
    if iFile < 5; %normal timestep
      y = y1; m = m1; d = d1; FilePath = FilePath1; Hour = (iFile-1)*6; 
    else %first of next day
      y = y2; m = m2; d = d2; FilePath = FilePath2; Hour = 0; 
    end

    
    File.T = [FilePath,'siganl.gdas.',sprintf('%04d',y),sprintf('%02d',m),...
              sprintf('%02d',d),sprintf('%02d',Hour),'.T.nc4'];
    File.p = [FilePath,'siganl.gdas.',sprintf('%04d',y),sprintf('%02d',m),...
              sprintf('%02d',d),sprintf('%02d',Hour),'.p.nc4'];
    Lat  = full(ncArray(File.T,'lat'));
    Lon  = full(ncArray(File.T,'lon'));
    T    = full(ncArray(File.T,'ta'));
    Time = datenum(y,m,d,Hour,0,0);
    p    = full(ncArray(File.p,'pfull'))./100;
    
    if iFile == 1;
      Data.latitude  = Lat;
      Data.longitude = Lon;
      Data.T         = T;
      Data.time      = Time;
      Data.p         = p;
    else
      Data.T    = cat(4,Data.T,T);
      Data.time = cat(1,Data.time,Time);
    end
    clear Lat Lon T Time p
  
  catch
    
    %failed to load a file. skip out
    Model.Error = 2;
    return
    
  end
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
PrsScale = squeeze(nanmean(nanmean(Data.p,1),2));

%then, stick stuff in a struct
Model.Lon  = Data.longitude;
Model.Lat  = Data.latitude;
Model.Time = Data.time;
Model.T    = permute(Data.T,[4,1,2,3]);
Model.Prs  = PrsScale;

%CFSR runs from -180 to +180, so we don't need to fix the lons

%success!
Model.Error = 0;
return


