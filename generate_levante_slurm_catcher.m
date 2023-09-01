clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%script to generate slurm scripts needed for DKRZ AIRS sampling runs
%
%Corwin Wright, c.wright@bath.ac.uk, 2023/04/06
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% set the instrument and model
Settings.Instrument = 'AIRS3D';
Settings.Model      = 'dyamond_icon5km';

%choose the data to sample
Settings.Dates   = datenum(2020,2,14);%:1:datenum(2020,3,1);
Settings.SubSets = 1:1:280;

%check if we have generated a track file for this, and only generate a script if we have
Settings.CheckExists = 1;
Settings.CheckDir = 'C:\Data\corwin\sampling_project\tracks\';

%how many jobs in each slurm call?
Settings.JobsPerCall = 20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make list of options, job names, and output file paths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[a,b] = ndgrid(Settings.Dates,Settings.SubSets);

Jobs.Dates   = a(:);
Jobs.SubSets = b(:);

Jobs.Names    = {};
for iJob=1:1:numel(a)
  Jobs.Names{   iJob} = [Settings.Instrument(1),Settings.Model(9:11),num2str(floor((a(iJob)/1000-floor(a(iJob)/1000))*1000))];  
end; clear a b c d e f iJob


%check file exists
if Settings.CheckExists == 1;

  for iJob=1:1:numel(Jobs.Names)
    FileToCheck = [Settings.CheckDir,'/',Settings.Instrument,'/track_',Settings.Instrument,'_', ...
                   num2str(Jobs.Dates(iJob)),'_g',sprintf('%03d',Jobs.SubSets(iJob)),'.mat'];
    if ~exist(FileToCheck,'file')
      %remove from the list
      Jobs.Names{iJob} = 'SKIP_THIS';
    end
  end

end

%sort by date
[~,idx] = sort(Jobs.Dates,'asc');
Jobs.Dates = Jobs.Dates(idx);
Jobs.SubSets = Jobs.SubSets(idx);
Jobs.Names = Jobs.Names(idx);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create a slurm script for each job
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


FireList = [];
for iJob = 1:Settings.JobsPerCall:numel(Jobs.Names)

  %check if we actually want this job
  if strcmp(Jobs.Names{iJob},'SKIP_THIS'); continue;
  else FireList = [FireList,iJob];
  end


  %create script header
  Script = "#!/bin/bash";
  
  Script(end+1) = "#SBATCH --account=bm1233";                %account name for budget.
  Script(end+1) = "#SBATCH --job-name="+Jobs.Names{iJob};    %name of job. Generated automatically.
  Script(end+1) = "#SBATCH --output=output.%j";              %standard (text) output file.
  Script(end+1) = "#SBATCH --error=error.%j";                %standard (text) error file.
  Script(end+1) = "#SBATCH --ntasks=16";                     %tasks allowed.     
  Script(end+1) = "#SBATCH --partition=compute";             %partition (queue) to use  
  Script(end+1) = "#SBATCH --time=8:00:00";                 %max allowed time

  %couple of lines of padding
  Script(end+1) = ""; Script(end+1) = ""; 

  %change to directory changing the code
  Script(end+1) = "cd ~/Code/Sampling";                      %wherever you put the code. CHANGE THIS.

  %couple of lines of padding
  Script(end+1) = ""; Script(end+1) = ""; 

  %load matlab and CDO modules
  Script(end+1) = "module load matlab";
  Script(end+1) = "module load cdo";

  %generate matlab commands
  Command = 'matlab -r "';
  for iCommand=iJob:1:iJob+Settings.JobsPerCall
    if iCommand > numel(Jobs.Names); continue; end
    Command = [Command,'[Error,OldData] = sampling_core_v2(''',Settings.Instrument,''',''',Settings.Model,''',',num2str(Jobs.Dates(iJob))];
    Command = [Command,',''Subset'',',num2str(Jobs.SubSets(iCommand)),',''TextUpdate'',true);'];
  end; clear iCommand

  %finalise the command and the script
  Command = [Command,'exit;"'];
  Script(end+1)  = Command;

  %write the script out as a text file
  writelines(Script,['job_',sprintf('%06d',iJob),'.txt'],LineEnding="\n");

  clear Sensitivity Command Script

end; clear iJob


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create a shell script to launch the jobs as a batch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Shell = "";
for iJob=FireList
  Shell(end+1) = "sbatch --dependency=singleton "+['job_',sprintf('%06d',iJob),'.txt'];
  Shell(end+1) = "rm "                           +['job_',sprintf('%06d',iJob),'.txt'];  
end
Shell(end+1) = "rm jobs_queue.sh";
writelines(Shell,'jobs_queue.sh',LineEnding="\n");
clear Shell iJob
