function CoreVars = sampling_core_variables

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  %common storage location for variables used in multiple files
  %note that variables which are physically-based are hardcoded - this file 
  %contains only variables which are nominally arbitrary, or at least in our 
  %control
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  %general settings
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %master file path for files generated
  CoreVars.MasterPath = [LocalDataDir,'/corwin/sampling_project/'];
 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  %model-specific settings
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %these have grown too numerous to be stored here
  %so they are now handled in the individual functions
  %just add the path where the functions are stored
  addpath('../models/')
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  %instrument-specific settings
  %most of these can be ovewritten in sensitivity-testing mode
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %COSMIC
  %%%%%%%%%%%%%%%%%%%%%%%%
  CoreVars.Cosmic.Path         = [LocalDataDir,'/COSMIC/daily_atmprf'];
  CoreVars.Cosmic.HeightRange  = h2p([8,50]);
  CoreVars.Cosmic.TimeRange    = [datenum(2006,1,112),datenum(2016,12,30)];  
  CoreVars.Cosmic.FineGrid.X   = 10; %along-track km. chosen from sensitivity tests
  CoreVars.Cosmic.FineGrid.Y   = 0.5; %across-track km. chosen from sensitivity tests
  CoreVars.Cosmic.FineGrid.Prs = 1/80; %decades of pressure. chosen from sensitivity tests
  
  %AIRS (2D - NASA JPL)
  %%%%%%%%%%%%%%%%%%%%%%%%
  CoreVars.Airs.Path         = [PhysicalDataDir,'/AIRS/selectedchannels/'];
  CoreVars.Airs.HeightRange  = [2.5]; %%[2,2.5,3,4,7,10,20,30,40,60,80,100]; 
  CoreVars.Airs.TimeRange    = [datenum(2002,8,31),datenum(2016,12,31)];    
  CoreVars.Airs.FineGrid.X   = 2; 
  CoreVars.Airs.FineGrid.Y   = 23;   
  CoreVars.Airs.FineGrid.Prs = 1/20;

  %AIRS (3D - Hoffman and Alexander)
  %%%%%%%%%%%%%%%%%%%%%%%%
  CoreVars.Airs3D.Path         = [LocalDataDir,'/AIRS/3d_airs/'];
  CoreVars.Airs3D.HeightRange  = h2p([15,65]); 
  CoreVars.Airs3D.TimeRange    = [datenum(2013,1,1),datenum(2013,12,31)];   
  CoreVars.Airs3D.FineGrid.X   = 1; 
  CoreVars.Airs3D.FineGrid.Y   = 1;  
  CoreVars.Airs3D.FineGrid.Prs = 1/20;
     
  %HIRDLS
  %%%%%%%%%%%%%%%%%%%%%%%%
  CoreVars.Hirdls.Path         = [LocalDataDir,'/HIRDLS/HIRDLS-L2'];
  CoreVars.Hirdls.HeightRange  = h2p([15,80]);
  CoreVars.Hirdls.TimeRange    = [datenum(2005,1,29),datenum(2008,3,31)];    
  CoreVars.Hirdls.FineGrid.X   = 10; 
  CoreVars.Hirdls.FineGrid.Y   = 2;   
  CoreVars.Hirdls.FineGrid.Prs = 1/80;
    
  %SABER
  %%%%%%%%%%%%%%%%%%%%%%%%
  CoreVars.Saber.Path         = [LocalDataDir,'/SABER/reprocessed-v2'];
  CoreVars.Saber.HeightRange  = h2p([15,80]);
  CoreVars.Saber.TimeRange    = [datenum(2002,1,1),datenum(2016,12,31)];
  CoreVars.Saber.FineGrid.X   = 10; 
  CoreVars.Saber.FineGrid.Y   = 2;  
  CoreVars.Saber.FineGrid.Prs = 1/80;
   
  
end

