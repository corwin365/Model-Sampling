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
  
  CoreVars.CESM_CK.Path   = [LocalDataDir,'/CESM/'];
  CoreVars.CESM_CK.Period = [datenum(2010,10,1),datenum(2010,10,31)]; 
  
  CoreVars.CFSR.Path   = [LocalDataDir,'/CFSR/'];
  CoreVars.CFSR.Period = [datenum(2002,1,1),datenum(2011,3,31)];
  
  CoreVars.ERAi.Path   = [LocalDataDir,'/ERAi/'];
  CoreVars.ERAi.Period = [datenum(2002,1,1),datenum(2012,12,31)];  
  CoreVars.ERAI = CoreVars.ERAi; %alias
  
  CoreVars.ERA5.Path   = [LocalDataDir,'/ERA5/'];
  CoreVars.ERA5.Period = [datenum(2010,1,1),datenum(2016,12,31)];  
  
  CoreVars.EC_FC.Path   = [LocalDataDir,'/ECMWF_fc'];
  CoreVars.EC_FC.Period = [datenum(2018,1,1),datenum(2018,1,365)];
  
  CoreVars.SONJA.Path   = [LocalDataDir,'/gisinger/'];
  CoreVars.SONJA.Period = [datenum(2018,1,92),datenum(2018,1,245)];  
  
  CoreVars.EcOpAl.Path   = [LocalDataDir,'/ECMWF/'];
  CoreVars.EcOpAl.Period = [datenum(2010,1,1),datenum(2011,1,31)];  
  
  CoreVars.JRA55.Path  = [LocalDataDir,'/JRA55/'];
  CoreVars.JRA55.Period = [datenum(2002,1,1),datenum(2012,12,31)];  
    
  CoreVars.JRA55C.Path = [LocalDataDir,'/JRA55C/'];
  CoreVars.JRA55C.Period = [datenum(2002,1,1),datenum(2012,12,31)];    
  
  CoreVars.MERRA2.Path = [LocalDataDir,'/MERRA2/'];
  CoreVars.MERRA2.Period = [datenum(2002,1,1),datenum(2016,12,31)];    
  
  CoreVars.UM_FC.Path = [LocalDataDir,'/corwin/annelize'];
  CoreVars.UM_FC.Period = [1,1].*datenum(2018,8,12);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  %instrument-specific settings
  %most of these can be ovewritten in sensitivity-testing mode
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %COSMIC
  %%%%%%%%%%%%%%%%%%%%%%%%
  CoreVars.Cosmic.Path         = [LocalDataDir,'/COSMIC/daily_atmprf'];
  CoreVars.Cosmic.HeightRange  = h2p([15,50]);
  CoreVars.Cosmic.TimeRange    = [datenum(2006,1,112),datenum(2016,12,30)];  
  CoreVars.Cosmic.FineGrid.X   = 10; %along-track km. chosen from sensitivity tests
  CoreVars.Cosmic.FineGrid.Y   = 0.5; %across-track km. chosen from sensitivity tests
  CoreVars.Cosmic.FineGrid.Prs = 1/80; %decades of pressure. chosen from sensitivity tests
  
  %AIRS
  %%%%%%%%%%%%%%%%%%%%%%%%
  CoreVars.Airs.Path         = [PhysicalDataDir,'/AIRS/selectedchannels/'];
  CoreVars.Airs.HeightRange  = [2.5]; %%[2,2.5,3,4,7,10,20,30,40,60,80,100]; 
  CoreVars.Airs.TimeRange    = [datenum(2002,8,31),datenum(2016,12,31)];    
  CoreVars.Airs.FineGrid.X   = 2; 
  CoreVars.Airs.FineGrid.Y   = 23;   
  CoreVars.Airs.FineGrid.Prs = 1/20;

  %AIRS 3D
  %%%%%%%%%%%%%%%%%%%%%%%%
  CoreVars.Airs3D.Path         = [LocalDataDir,'/AIRS/3d_airs/'];
  CoreVars.Airs3D.HeightRange  = h2p([15,65]); 
  CoreVars.Airs3D.TimeRange    = [datenum(2013,1,1),datenum(2013,12,31)]; 
% % %   %not yet sensitivity-tested    
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

