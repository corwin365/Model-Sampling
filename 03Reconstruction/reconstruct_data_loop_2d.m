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
Instruments = {'AIRS_1D','AIRS_QBO','AIRS'};
Models      = {'ERA5','CFSR','JRA55','JRA55C','ERAI','MERRA2'};

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
      
      InFile = FileList{iFile};
      
      OutFile = strfind(InFile, 'sampled');
      OutFile = [Settings.OutDataDir,'/',Settings.Instrument,'/',Settings.Model,'/',InFile(OutFile:end)];
      if exist(OutFile,'file') ~= 0;
%         disp([OutFile,' already done']);
       continue
      end
      
      disp([OutFile,' started']);
      
      %% load the data
      Data = load(InFile);
      
      if isfield(Data.Output,'T');
        
        
        %% preprocess
        %first, work out the maximum points per day
        Ng =  max(Data.Recon.g);
        Nx  = max(Data.Recon.x);
        Ny  = max(Data.Recon.y);
        %and the maximum number of levels
        NLevels = max(Data.Recon.z);

        %hence, create arrays to store the reconstructed data
        Sampled_Data         = struct();
        Sampled_Data.T       = NaN(Ng,Nx,Ny,NLevels);
        Sampled_Data.Tsimple = Sampled_Data.T;  
        Sampled_Data.Lat     = Sampled_Data.T;
        Sampled_Data.Lon     = Sampled_Data.T;
        Sampled_Data.Time    = Sampled_Data.T;
        Sampled_Data.Prs     = Sampled_Data.T;
        
        %% extract
        taxis = floor(Data.Output.Time-floor(min(Data.Output.Time(:))))+1;
        
        
        for iPoint=1:1:numel(Data.Output.Lat);
          
          g = Data.Recon.g(iPoint);
          x = Data.Recon.x(iPoint);
          y = Data.Recon.y(iPoint);      
          z = Data.Recon.z(iPoint);
          t = taxis(iPoint);
          
          Sampled_Data.T(      g,x,y,z) = Data.Output.T(   iPoint);
          Sampled_Data.Tsimple(g,x,y,z) = Data.Output.TSimple(   iPoint);          
          Sampled_Data.Lat(    g,x,y,z) = Data.Output.Lat( iPoint);
          Sampled_Data.Lon(    g,x,y,z) = Data.Output.Lon( iPoint);
          Sampled_Data.Time(   g,x,y,z) = Data.Output.Time(iPoint);
          Sampled_Data.Prs(    g,x,y,z) = 10.^Data.Output.Prs( iPoint); %we worked in logs
          
        end
        
        %we already have a time axis, so no need to split by day
        Vars = {'T','Lat','Lon','Time','Prs'};


        %rename and save
        parsave(OutFile,Sampled_Data);
        disp([OutFile,' done']);
        
      else
        %input file: skip
        disp([OutFile,' skipped: no temperature data']);
      end
      
    end
    
  end
end
