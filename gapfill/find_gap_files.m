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

%paths
Settings.Dir.Output = [LocalDataDir,'/corwin/sampling_project/output/AIRS3D/era5_merged/test'];
Settings.Dir.Era5   = [LocalDataDir,'/ERA5/fullfat/'];
Settings.Dir.Airs   = [LocalDataDir,'/AIRS/3d_airs'];

%dates to consider
Settings.TimeScale =  datenum(2007,1,1:9);%1:1:365);%:1:datenum(2020,2,1);

%output filenames
Settings.OutFile.Airs = 'airs_missing.mat';
Settings.OutFile.Era5 = 'get_era5.py';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create storage arrays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ToGet.Airs = [];
ToGet.Era5 = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% process!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iDay=1:1:numel(Settings.TimeScale)

  %% find missing granules on this day
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  MissingList = zeros(240,1);

  %load the storage file for this day
  StoreFile = [Settings.Dir.Output,'/',yyyyDDD(Settings.TimeScale(iDay)),'.mat'];
  if ~exist(StoreFile,'file')
    %all granules are 'missing' - add to list for this day
    MissingList(:) = 1;
  else
    %load the file and identify the missing granules
    try
      StoreFile = load(StoreFile,'T');
      MissingList(find(nansum(StoreFile.T,[2,3,4])  == 0)) = 1;
    catch; 
      MissingList(:) = 1;
    end
  end
  clear StoreFile

  %if we have no missing granules, we're happy
  if sum(MissingList) == 0; clear MissingList; 
    disp([datestr(Settings.TimeScale(iDay)),' fully processed'])
    continue; 
  end

  %% see if we have AIRS data for the missing files
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %get a list of missing files
  MissingList = find(MissingList == 1);


  %try and load them all
  ShouldBeThere = zeros(size(MissingList));
  for iG=1:1:numel(MissingList)
    [y,~,~] = datevec( Settings.TimeScale(iDay));
    dn      = date2doy(Settings.TimeScale(iDay));
    AirsFile = [Settings.Dir.Airs,'/',sprintf('%04d',y),'/',sprintf('%03d',dn),'/', ...
               'airs_',sprintf('%04d',y),'_',sprintf('%03d',dn),'_',sprintf('%03d',MissingList(iG)),'.nc'];
    if exist(AirsFile,'file'); ShouldBeThere(iG) = 1; end
  end
  clear iG y dn AirsFile

  %if all the missing granules from our data are also missing from AIRS storage, we're stillhappy
  if sum(ShouldBeThere) == 0;
    clear MissingList ShouldBeThere;
    disp([datestr(Settings.TimeScale(iDay)),' fully processed'])
    continue; 
  end

  %% ok, the granule is actually missing. What do we need
  % to download to fix this
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  MissingList = MissingList(ShouldBeThere == 1);
  clear ShouldBeThere

  %we need a list of missing granules. This is easy
  ToGet.Airs = cat(1,ToGet.Airs,[ones(size(MissingList)).*Settings.TimeScale(iDay),MissingList]);

  %we also need the ERA5 hourly timesteps surrounding them. This is fiddlier.
  %first, get a list of the hours we need
  BackgroundList = NaN(numel(MissingList).*2,3);
  [y,~,~] = datevec( Settings.TimeScale(iDay));
  dn      = date2doy(Settings.TimeScale(iDay));
  for iG=1:1:numel(MissingList)

    %find the hours during and after the granule
    hh = floor(MissingList(iG)./10);
    BackgroundList(2.*(iG-1)+1,:) = [y,dn,hh];
    BackgroundList(2.*(iG-1)+2,:) = [y,dn,hh+1];
  end; clear y dn hh iG MissingList

  %adjust entries that fall on the day after
  DayAfter  = find(BackgroundList(:,3) >= 24);
  BackgroundList(DayAfter,2) = BackgroundList(DayAfter,2)+1;
  BackgroundList(DayAfter,3) = BackgroundList(DayAfter,3)-24;
  clear DayAfter

  %adjust entries that fall on the day before
  DayBefore  = find(BackgroundList(:,3) < 0);
  BackgroundList(DayBefore,2) = BackgroundList(DayBefore,2)-1;
  BackgroundList(DayBefore,3) = BackgroundList(DayBefore,3)+24;
  clear DayBefore

  %drop this down to a list of unique entries
  BackgroundList = unique(BackgroundList,'rows');

  %and store
  
  ToGet.Era5 = [ToGet.Era5;BackgroundList];
  disp([datestr(Settings.TimeScale(iDay)),' missing ',num2str(numel(find(ToGet.Airs(:,1) == Settings.TimeScale(iDay)))),' granules; added to list.'])

  clear BackgroundList

end; clear iDay

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%do we have any missing granules?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if numel(ToGet.Era5) == 0

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %we're done
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  disp('No missing granules found!')

else

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% write the results of this to two files:
  %
  % 2. a Matlab save file for the AIRS granule list
  % 1. a Python download file for the ERA5 data
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %Python download script
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %generate needed variables
  YearsOut = "YearList = [";
  DaysOut  = "DaysList = [";
  HoursOut = "HourList = [";
  for iFile=1:1:size(ToGet.Era5,1)
    if iFile == size(ToGet.Era5,1); EndChar = "]"; else EndChar = ","; end
    YearsOut = YearsOut+num2str(ToGet.Era5(iFile,1))+EndChar;
    DaysOut  =  DaysOut+num2str(ToGet.Era5(iFile,2))+EndChar;
    HoursOut = HoursOut+num2str(ToGet.Era5(iFile,3))+EndChar;
  end
  clear iFile EndChar

  %put them in a script
  Script        = "#!/usr/bin/env python";
  Script(end+1) = "import os";
  Script(end+1) = "import sys";
  Script(end+1) = "import datetime";
  Script(end+1) = "import numpy";
  Script(end+1) = "import math";
  Script(end+1) = "import cdsapi";
  Script(end+1) = "from joblib import Parallel, delayed";
  Script(end+1) = "from time import sleep";
  Script(end+1) = "from random import randint";
  Script(end+1) = "";
  Script(end+1) = "";
  Script(end+1) = "#disable unnecessary warnings";
  Script(end+1) = "import urllib3";
  Script(end+1) = "urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)";
  Script(end+1) = "";
  Script(end+1) = "";
  Script(end+1) = "#declare function to get data for each day";
  Script(end+1) = "def getECMWF(Year,Day,Hour):";
  Script(end+1) = "  ";
  Script(end+1) = "  Year         = int(Year)";
  Script(end+1) = "  Day          = int(Day)";
  Script(end+1) = "  Hour         = int(Hour)";
  Script(end+1) = "";
  Script(end+1) = "";
  Script(end+1) = "  FilePath = '"+strrep(Settings.Dir.Era5,'\','\\')+"'";
  Script(end+1) = "  FileName = 'era5_' + format(Year, '04d') + 'd' + format(Day,  '03d') + 'h' + format(Hour,'02d')";
  Script(end+1) = "  FileName += '.nc'";
  Script(end+1) = "";
  Script(end+1) = "  DateString = format(Year, '04d')+'-'+format(1, '02d')+'-'+format(Day, '02d')";
  Script(end+1) = "  ";
  Script(end+1) = "  DateString = datetime.datetime(Year,1,1)+datetime.timedelta(days=Day-1)";
  Script(end+1) = "  DateString = DateString.strftime('%Y-%m-%d')";
  Script(end+1) = "";
  Script(end+1) = "  print('================================================')";
  Script(end+1) = "  print(FilePath+FileName)";
  Script(end+1) = "  print('================================================')";
  Script(end+1) = "";
  Script(end+1) = "";
  Script(end+1) = "    ";
  Script(end+1) = "  if not os.path.exists(FilePath+FileName):";
  Script(end+1) = "";
  Script(end+1) = "      ##space out requests";
  Script(end+1) = "      #sleep(randint(10,180))";
  Script(end+1) = "";
  Script(end+1) = "      #try:";
  Script(end+1) = "        #get the data!";
  Script(end+1) = "      c = cdsapi.Client()";
  Script(end+1) = "      c.retrieve('reanalysis-era5-complete', {    # do not change this!";
  Script(end+1) = "                'class': 'ea',";
  Script(end+1) = "                'dataset': 'era5',";
  Script(end+1) = "                'date': DateString,";
  Script(end+1) = "                'expver': '1',";
  Script(end+1) = "                'levelist': '1/to/90',";
  Script(end+1) = "                'levtype': 'ml',";
  Script(end+1) = "                'param': '130',";
  Script(end+1) = "                'stream': 'oper',";
  Script(end+1) = "                'time':  format(Hour,'02d')+':00:00',";
  Script(end+1) = "                'type': 'an',";
  Script(end+1) = "                'format': 'netcdf',";
  Script(end+1) = "                'grid': '0.25/0.25',";
  Script(end+1) = "                'area': '90/0/-90/360',";
  Script(end+1) = "                }, FilePath+FileName)";
  Script(end+1) = "        ";
  Script(end+1) = "      #except:";
  Script(end+1) = "        #print('Error downloading file '+FileName)";
  Script(end+1) = "      ";
  Script(end+1) = "  return";

  Script(end+1) = YearsOut;
  Script(end+1) = DaysOut;
 Script(end+1) = HoursOut; 

 Script(end+1) = "results = Parallel(n_jobs=10, verbose=11, backend='threading')(map(delayed(getECMWF),YearList,DaysList,HourList))";


 %output
 fileID = fopen(Settings.OutFile.Era5,'w');
 for iLine=1:1:numel(Script)
   Line = strrep(Script{iLine},'%','%%');
   fprintf(fileID, Line);
   fprintf(fileID, '\n');   
 end
 fclose(fileID);
 clear fileID iLine Script Line

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %get ERA5 and produce trackfiles
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %produce track files
 disp([num2str(size(ToGet.Airs,1))+" missing AIRS granules; generating"])
 for iG=1:1:size(ToGet.Airs,1)
   oif_airs_function(ToGet.Airs(iG,1),ToGet.Airs(iG,2))
 end

 %also write out the list to a sae file so we can sample it
 Needed = ToGet.Airs;
 save(Settings.OutFile.Airs,'Needed')
 disp(["List of "+num2str(size(Needed,1))+" missing AIRS granules written to file"])
 clear Needed


 %fire off Python script to get ERA5 data
 disp(['Using Python to download ',num2str(size(ToGet.Era5,1)),' ERA5 hourly files'])
 system("python3 "+Settings.OutFile.Era5)

end
