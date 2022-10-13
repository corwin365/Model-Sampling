clearvars
addpath(genpath('../common'));
CoreVars = sampling_core_variables;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%reconstruct data extracted from model
%works for 3D tracks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% settings

%input directory
Settings.InDataDir = [CoreVars.MasterPath,'/samples/'];

%instrument
Instruments = {'AIRS3D'};
Models      = {'IFS_1KM_FAKEWAVE'}; %


for iInst=1:1:numel(Instruments)
  for iModel=1:1:numel(Models);
   disp([Instruments{iInst},'   ',Models{iModel}])
    
    %instrument
    Settings.Instrument = Instruments{iInst};
    
    %model
    Settings.Model = Models{iModel};
    
    %hence, path to files
    Settings.DataDir = [Settings.InDataDir,'/',Settings.Instrument,'/',Settings.Model];
    
    %get a list of all the tracks in that directory
    warning off
    FileList = wildcardsearch(Settings.DataDir,'*.mat');
    warning on
    
    %loop over days
    for iFile=1:1:numel(FileList);
           
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %load the data
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      %output filename. 
      OutFile = strrep(FileList{iFile},'samples','reconstructed');

      %check it doesn't already exist
      if exist(OutFile,'file');
%          disp([OutFile,' already done']);
        continue
      end

      Data = load(FileList{iFile});

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %ok, regrid!
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      Data.Recon.g = ones(size(Data.Recon.x));
      
      %reorder the data into the desired output sequence and then reshape.
      Nx = max(Data.Recon.x);
      Ny = max(Data.Recon.y);
      Nz = max(Data.Recon.z);
      Ng = max(Data.Recon.g);
      Variables = {'T','TSimple','Lat','Lon','Time','Prs'};
      [~,Order]  = sortrows([Data.Recon.g,Data.Recon.z,Data.Recon.y,Data.Recon.x],[1,2,3,4]);

      Sampled_Data = struct();
      
      for iVar=1:1:numel(Variables);
        
        Var = Data.Output.(Variables{iVar});
        Var = Var(Order);
        Var = reshape(Var,Nx,Ny,Nz,Ng);
        
        Sampled_Data.(Variables{iVar}) = Var;
      end
        
        
      parsave(OutFile,Sampled_Data);
      disp([OutFile,' done']);
    end
  end
end
