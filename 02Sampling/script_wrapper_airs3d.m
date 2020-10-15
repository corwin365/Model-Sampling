clearvars
addpath('../common');
CoreVars = sampling_core_variables;

%parameters to vary
ToVary.Instruments = {'airs3d'}%instruments to process
ToVary.Models      = {'CFSR','JRA55','MERRA2'}; %models to process

%node runtime for a given instrument, minutes. programme wil select this.
%getting this right allows more parallel jobs to run inside balena group cputime limit
RunTime.airs3d = 360;
RunTime.gisinger = 60;
RunTime.airs_qbo = 360;

%how many granules should we try and process in each pass?
GranuleStep = 80; %number of granules per pass

%forecast? if so, what time step?
%set to 0 to use non-forecast data
ForecastHours = 72;

%clobber?
NoClobber  = 0;

Count = 1;
for iModel=1:1:numel(ToVary.Models);
  for iInst=1:1:numel(ToVary.Instruments);
    for iYear=2018:1:2018
      for iDay=224:230
        for iGranule=1:GranuleStep:240
        
          GranRange = [iGranule,iGranule+GranuleStep-1];

%            Ins = ToVary.Instruments{iInst}; Ins = Ins(1:2);
          Mod = ToVary.Models{iModel}; Mod = [upper(Mod(1)),upper(Mod(end))];
          YY  = sprintf('%02d',iYear-2000);
          DD  = sprintf('%03d',iDay);
          
          JobName = [Mod,YY,DD,sprintf('%03d',iGranule)];
          


          Time = num2str(RunTime.(ToVary.Instruments{iInst}));
          
            
          %ERA5 special handling
          if strcmp(ToVary.Models{iModel},'era5') == 1 | strcmp(ToVary.Models{iModel},'ec_fc');          
            if iYear < 2010; continue; end
            NThreads  = 4;
            Partition = 'batch-sky'; %only the Skylake nodes are ready for this jelly
            Time = num2str(360);
            CPUSPerNode = 24;
            
          else
            Partition = 'batch-short';
            CPUSPerNode = 16;  
            NThreads = 16;
          end
          
          
          Text1 = ['#!/bin/bash\n## Name of the job\n#SBATCH --job-name=',JobName,'\n## Account to charge to\n#SBATCH --account=free\n\n'];
          Text2 = ['\n#SBATCH --time=',Time,':00\n## Number of node required and tasks per node\n'];
          Text3 = ['#SBATCH --nodes=1\n#SBATCH --ntasks-per-node=',num2str(CPUSPerNode),'\n\n#SBATCH --output=log.%%j.out\n#SBATCH --error=log.%%j.err\n#SBATCH --partition=',Partition];
          
          
          

        

                  
          
          
          %Load Matlab environment
          Text4 = ['\n\n\nmodule load matlab'];
          
          
          Text5 = ['\nmatlab -nodisplay -r "'];
          
          
          %i'm getting a problem with parallel pools not loading
          %I think it might be due to Balena trying too create many at the same time
          %so add a random wait of between 0-30 seconds to each script
          Wait = rand.*30;
          Commands = ['pause(',num2str(Wait),');'];
          
            
          Commands = [Commands,'Instruments = {''',ToVary.Instruments{iInst},'''};TimeScale = [datenum(',num2str(iYear),',1,',num2str(iDay),'):1:datenum(',num2str(iYear),',1,',num2str(iDay),')];NoClobber=',num2str(NoClobber),';Granules=[',num2str(GranRange(1)),',',num2str(GranRange(2)),'];ForecastHours=',num2str(ForecastHours),';Models = {''',ToVary.Models{iModel},'''};parpool(',num2str(NThreads),');wrapper'];

          Text6 = [';exit"'];
          
          fid = fopen(['job',sprintf('%04d',Count),'_wrapper.txt'],'wt');
          fprintf(fid, Text1);
          fprintf(fid, Text2);
          fprintf(fid, Text3);
          fprintf(fid, Text4);
          fprintf(fid, Text5);
          fprintf(fid, Commands);
          fprintf(fid, Text6);
          fclose(fid);
          
          Count = Count+1;
        end
      end
    end
  end
end


%generate file to fire the whole lot off
fid = fopen(['fire_wrappers.sh'],'wt');
%  for i=1:1:Count-1;fprintf(fid,['sbatch --begin=now+7hour job',sprintf('%04d',i),'.txt\n']);end
for i=1:1:Count-1;
   fprintf(fid,['sbatch  job',sprintf('%04d',i),'_wrapper.txt\n']);
%    fprintf(fid,['sbatch --dependency singleton job',sprintf('%04d',i),'_wrapper.txt\n']);  
%    fprintf(fid,['sbatch --dependency singleton job',sprintf('%04d',i),'_wrapper.txt\n']);    
  fprintf(fid,['rm job',sprintf('%04d',i),'_wrapper.txt\n']);    
end
fprintf(fid,['rm fire_wrappers.sh\n']); 
fclose(fid);

disp(['Written ',num2str(Count),' files (probably)'])

