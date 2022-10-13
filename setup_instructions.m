clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%manual setup instructions are in the first part of this script
%once you've followed them, run this script to create the necessary paths
%
%Corwin Wright, c.wright@bath.ac.uk, 2022/10/13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Instructions. READ THESE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%understanding the context
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%For an explanation of what's going on see Wright and Hindley, ACP 2018, doi: 10.5194/acp-18-13703-2018
%
% in terms of Figure A1 of that paper:
%  - the contents of './models/' are the "Model Import Functions" (MIFs)
%  - the contents of './01CreateTracks/' are the 'Observation Import Functions' (OIFs)
%  - the contents of './02Sampling/' are the "Core Analysis"
%  - the contents of './03Reconstruction/' are used to put the reformat the output to be human-readable
%  - the contents of './common/' are functions used in the above 
%
%
%It should probably go without saying that you need to be using input data (model and obs) in exactly the same
%format I used when I wrote the relevant functions. If you don't have this, either ask me to tell you what it was, or write
%your own OIFs/MIFs as needed. As long as they have consistent inputs and outputs with the existing ones they should work fine, I 
%have added several over the years without needing to modify the base code
%
%  - for OIFs, just stick them in the 01CreateTracks directory and add the directory they OUTPUT to the list in ./common/oif_list.m
%  - for MIFs, stick them in ./models, but you also need to add their name and list of arguments to ./common/model_settings.m 
%    so they can be called from the main sampling routine when needed.


%What you need to do to get it working on your system:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%note: in all paths of my base version, the function 'LocalDataDir' is just a function I use locally to help manage my filesystem
%You shouldn't include this in your manually set paths for 1a, 1b, 2a and 3a - absolute paths from root are probably easiest

%1. in ./common/sampling_core_variables.m:
%  a. set CoreVars.MasterPath to the root of the directory where you want to store output
%  b. set CoreVars.[Instrument Name].Path for each instrument you want to use to where you have raw input data from that instrument stored
%
%
%2. For each model you wish to sample from, find the relevant function in ./models and:
%  a. set CoreVars.[Model Name].Path to where you have the relevant model fields stored  on your file system
%(if you're not sure what the function you want is, there's a listing in ./common/model_settings saying which model needs which function)
%
%
%3. For each observation type you want to sample as, find the relevant function in ./01CreateTracks and:
%  a. set Settings.InDir to the location of the raw observational data on your file system
%
%
%4. Once that's done, run this programme and it should create all the necessary subdirectories in CoreVars.MasterPath

%What you need to do to run it:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%1. use the functions in ./01CreateTracks to create track-parameter files for your chosen instrument and period
%2. in ./02Sampling, use 'wrapper.m' or 'wrapper_singlethread.m' to process the data
%3. once the data are processed, use the appropriate file in ./03Reconstruction to make the data format more readable (optional)
%
%for all of them, you will need to set various options - these options are always at the top of the programme, I never bury them further down
%


%for (2), the only difference between the functions is that one uses Matlab's inbuild parallelisation libraries to run
%loads of samples at the same time from the same call in 'wrapper.m', while the singlethread version only does one at a time
% the standard version will be faster for most datasets, but for high-resolution datasets the singlethread version uses less 
%memory (the multi version uses roughly NThreads x as much), and hence is more likely to run successfully
%

%for 3, I use:
% - reconstruct_data_loop for 1D datasets like limb sounders or radionsondes
% - reconstruct_data_loop_2d for 2D datasets like single AIRS channels or airglow imagery
% - reconstruct_data_loop_3d_v2 for 3D datasets like AIRS, when procvessing entire days globally
% - reconstruct_data_loop_3d_v2_withgaps for 3D datasets like AIRS, when the data has gaps in relative to the daily total (for example because we're only working regionally)
%
%You do not necessarly need to do this - the output from step (3) is perfectly correct but is jsut a lsit of points rather than formatted like the instrument measurement tracks


%if you only do step (2), then output will apear in the appropriate subdirectory of samples/ 
%if you do step (3) as well, then the above is true but there will also be reformatted data in reconstructed/



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Once you've FOLLOWED the instructions, run this file as a programme
%and the below will set up the system paths you need
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%get master path 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%add directory containing common functions to path.
warning('off','MATLAB:mpath:nameNonexistentOrNotADirectory')
addpath('./common')
warning('on','MATLAB:mpath:nameNonexistentOrNotADirectory')

%call core variable file to get the path set above
CV = sampling_core_variables;
MasterPath = CV.MasterPath;
clear CV

%%get list of observational datasets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get list of obs datasets
ObsList = oif_list;


%%get model list
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Models = model_settings([],'JUSTLISTING');


%%hence, create all the necessary directories!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%master path
if exist(MasterPath,'dir'); disp([MasterPath,' already exists, skipping']);
else; mkdir(MasterPath);    disp(['Created ',MasterPath]);
end

%subdirectory paths
SubDirs.Obs  = {'tracks'};                  %paths that only have subdirs named after obs
SubDirs.Both = {'samples','reconstructed'}; %paths that need both model and obs info

%obs-only
for iDir = 1:1:numel(SubDirs.Obs)
  for iObs=1:1:numel(ObsList)
    Path = [MasterPath,'/',SubDirs.Obs{iDir},'/',ObsList{iObs},'/'];
    if exist(Path,'dir');   disp([Path,' already exists, skipping']);
    else;      mkdir(Path); disp(['Created ',Path]);
    end
  end
end

%obs-and-model
for iDir = 1:1:numel(SubDirs.Both)
  for iObs=1:1:numel(ObsList)
    for iMod=1:1:numel(Models);
      Path = [MasterPath,'/',SubDirs.Both{iDir},'/',ObsList{iObs},'/',upper(Models{iMod}),'/'];
      if exist(Path,'dir');   disp([Path,' already exists, skipping']);
      else;      mkdir(Path); disp(['Created ',Path]);
      end
    end
  end
end
