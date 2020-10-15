clearvars
addpath(genpath('../common'));
CoreVars = sampling_core_variables;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%reconstruct data extracted from model
%works for 2D tracks, e.g. HIRDLS and SABER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% settings

%input directory
Settings.InDataDir = [CoreVars.MasterPath,'/samples/'];

%out-path, generated from in-path
Settings.OutDataDir = [CoreVars.MasterPath,'/reconstructed/'];

%instrument
Instruments = {'SABER'};%{'COSMIC','HIRDLS','SABER'};
Models      = {'EC_FC'};%{'CFSR','JRA55','JRA55C','ERAI','MERRA2','ERA5'};

for iInst=1:1:numel(Instruments)
for iModel=1:1:numel(Models);  
  
Settings.Instrument = Instruments{iInst};

%model
Settings.Model = Models{iModel};


%hence, path to files
Settings.DataDir = [Settings.InDataDir,'/',Settings.Instrument,'/',Settings.Model];


%get a list of all the files we want to run it over
warning off
FileList = wildcardsearch(Settings.DataDir,'*.mat');
warning on

%loop over files
parfor iFile=1:1:numel(FileList);
  try
  InFile = FileList{iFile};
  
  OutFile = strfind(InFile, 'sampled_');
  OutFile = [Settings.OutDataDir,'/',Settings.Instrument,'/',Settings.Model,'/',InFile(OutFile:end)];
  if exist(OutFile,'file') ~= 0; 
%     disp([OutFile,' already done']);
    continue
  end
  
%    disp([OutFile,' started']);

  %% load the data
  Data = load(InFile);

  if isfield(Data.Output,'T'); 

    
    %% preprocess
    %first, work out the maximum profiles per day
    NProfs  = max(Data.Recon.x);
    %and the maximum number of levels
    NLevels = max(Data.Recon.z);

    
    %hence, create arrays to store the reconstructed data
    Sampled_Data         = struct();
    Sampled_Data.T       = NaN(NProfs,NLevels);
    Sampled_Data.Tsimple = Sampled_Data.T;    
    Sampled_Data.Lat     = Sampled_Data.T;
    Sampled_Data.Lon     = Sampled_Data.T;
    Sampled_Data.Time    = Sampled_Data.T;
    Sampled_Data.Prs     = Sampled_Data.T;
    
    %% extract
    taxis = floor(Data.Output.Time-floor(min(Data.Output.Time(:))))+1;
    
    
    for iPoint=1:1:numel(Data.Output.Lat);
      
      x = Data.Recon.x(iPoint);
      z = Data.Recon.z(iPoint);
      t = taxis(iPoint);
      
      Sampled_Data.T(      x,z) = Data.Output.T(      iPoint);
      Sampled_Data.Tsimple(x,z) = Data.Output.TSimple(iPoint);      
      Sampled_Data.Lat(    x,z) = Data.Output.Lat(    iPoint);
      Sampled_Data.Lon(    x,z) = Data.Output.Lon(    iPoint);
      Sampled_Data.Time(   x,z) = Data.Output.Time(   iPoint);
      Sampled_Data.Prs(    x,z) = 10.^Data.Output.Prs(iPoint); %we worked in logs
      
    end
    
    %we already have a time axis, so no need to split by day
    %we also may as well drop empty profiles caused by this split
    Vars = {'T','Lat','Lon','Time','Prs'};
   
    for iVar=1:1:numel(Vars); 
      sz = size(Sampled_Data.(Vars{iVar}));
      Sampled_Data.(Vars{iVar}) = reshape(Sampled_Data.(Vars{iVar}),sz(1),sz(2));
    end

    %rename and save
    parsave(OutFile,Sampled_Data);
    disp([OutFile,' done']);
  
  else
    %input file: skip
    disp([OutFile,' skipped: no temperature data']);
  end
  

catch;end
end
end
end
