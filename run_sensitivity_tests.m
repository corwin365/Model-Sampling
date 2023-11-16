clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%script to generate slurm scripts needed for sensitivity 
%testing of the sampling core for a given instrument
%
%Corwin Wright, c.wright@bath.ac.uk, 2023/03/21
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% set the instrument and model
Settings.Instrument = 'AIRS3D';
Settings.Model      = 'dyamond_geos3km';

%choose the data to use as a testbed
Settings.Date = datenum(2020,1,23);
Settings.SubSet = 201%3; %e.g. AIRS granules. Set to NaN if this is not relevant.

%% set the output path
Settings.OutRoot = '/sens/';

%% set the range of values to vary over




%grid values

%merged limb sounders
% Settings.FineGrid.X   = 1:1:10;    %averaging blob x-dimension, km
% Settings.FineGrid.Y   = 1:1:10;    %averaging blob y-dimension, km
% Settings.FineGrid.Prs = 1./[80:-10:20];  %averaging blob z-dimension, decades of pressure

%AIRS3D
Settings.FineGrid.X   = fliplr(logspace(log10(0.25),log10(10),7));    %averaging blob x-dimension, km
Settings.FineGrid.Y   = fliplr(logspace(log10(0.25),log10(10),7));    %averaging blob y-dimension, km
Settings.FineGrid.Prs = 1./fliplr(linspace(10,40,7));  %averaging blob z-dimension, decades of pressure

% % %blob values
% % Settings.BlobScale = 3; %number of standard deviations to compute sensitivity out to (+- from centre)
% % Settings.MinSignal = 0.99; %fraction of signal to require computation over
% % 
% % %lowest altitude of model data used
% % Settings.MaxPrs = 1000;


%time permitted
Settings.Time = "4:00:00";

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make list of options, job names, and output file paths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[a,b,c,d,e,f] = ndgrid(Settings.FineGrid.X,Settings.FineGrid.Y,Settings.FineGrid.Prs, ...
                       Settings.BlobScale,Settings.MinSignal,Settings.MaxPrs);

Jobs.FineGrid.X   = a(:);
Jobs.FineGrid.Y   = b(:);
Jobs.FineGrid.Prs = c(:);
Jobs.BlobScale    = d(:);
Jobs.MinSignal    = e(:);
Jobs.MaxPrs       = f(:);

Jobs.OutPaths = {};
Jobs.Names    = {};
for iJob=1:1:numel(a)
  Jobs.OutPaths{iJob} = [Settings.OutRoot,'/senstest_', ...
                         Settings.Instrument,'_',Settings.Model, ...
                         '_test',sprintf('%06d',iJob),'.mat'];
  Jobs.Names{   iJob} = ['sens',sprintf('%06d',iJob),'.mat'];  
end; clear a b c d e f iJob

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create a slurm script for each job
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iJob = 1:1:numel(Jobs.OutPaths)


  %create script header
  Script = "#!/bin/bash";
  
  Script(end+1) = "#SBATCH --account=xxxxxx";                %account name for budget.
  Script(end+1) = "#SBATCH --job-name="+[Settings.Instrument(1),Settings.Model(9:11),num2str(iJob)]; %name of job. Generated automatically.
  Script(end+1) = "#SBATCH --output=output.%j";              %standard (text) output file.
  Script(end+1) = "#SBATCH --error=error.%j";                %standard (text) error file.
  Script(end+1) = "#SBATCH --ntasks=16";                     %tasks allowed.     
  Script(end+1) = "#SBATCH --partition=compute";             %partition (queue) to use  
  Script(end+1) = "#SBATCH --time="+Settings.Time;           %max allowed time
  % Script(end+1) = "#SBATCH --mem=480GB";                     %memory


  %couple of lines of padding
  Script(end+1) = ""; Script(end+1) = ""; 

  %change to directory changing the code
  Script(end+1) = "cd ~/Code/Sampling";                      %wherever you put the code. CHANGE THIS.

  %couple of lines of padding
  Script(end+1) = ""; Script(end+1) = ""; 

  %load matlab module
  Script(end+1) = "module load matlab";
  Script(end+1) = "module load cdo";

  %generate sensitivity testing struct in the file
  Sensitivity = 'Sensitivity = struct(); ';
  Sensitivity = [Sensitivity,  'Sensitivity.FineGrid.X = ',num2str(Jobs.FineGrid.X(  iJob)),';'];
  Sensitivity = [Sensitivity,  'Sensitivity.FineGrid.Y = ',num2str(Jobs.FineGrid.Y(  iJob)),';'];
  Sensitivity = [Sensitivity,'Sensitivity.FineGrid.Prs = ',num2str(Jobs.FineGrid.Prs(iJob)),';'];  
  % Sensitivity = [Sensitivity,   'Sensitivity.BlobScale = ',num2str(Jobs.BlobScale(   iJob)),';'];  
  % Sensitivity = [Sensitivity,   'Sensitivity.MinSignal = ',num2str(Jobs.MinSignal(   iJob)),';'];  
  % Sensitivity = [Sensitivity,      'Sensitivity.MaxPrs = ',num2str(Jobs.MaxPrs(      iJob)),';'];  

  %path to write output
  Sensitivity = [Sensitivity,'Sensitivity.NewPath = ''',Jobs.OutPaths{iJob}, ''';'];  

  %generate matlab base command
  Command = ['matlab -r "',Sensitivity,'[Error,OldData] = sampling_core_v3(''',Settings.Instrument,''',''',Settings.Model,''',',num2str(Settings.Date)];
  if ~isnan(Settings.SubSet); Command =[Command,',''Subset'',',num2str(Settings.SubSet)]; end
  Command = [Command,', ''Sensitivity'',Sensitivity'];

  %finalise the command and the script
  Command = [Command,',''SensTestMode'',true);exit;"'];
  Script(end+1)  = Command;

  %write the script out as a text file
  writelines(Script,['job_',sprintf('%06d',iJob),'.txt'],LineEnding="\n");

  clear Sensitivity Command Script
 
end; clear iJob

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create a shell script to launch the jobs as a batch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Shell = "sbatch "+['job_',sprintf('%06d',1),'.txt'];
for iJob=2:1:numel(Jobs.Names); 
  Shell(end+1) = "sbatch "+['job_',sprintf('%06d',iJob),'.txt'];
end
writelines(Shell,'fire_jobs.sh',LineEnding="\n");
clear Shell iJob
