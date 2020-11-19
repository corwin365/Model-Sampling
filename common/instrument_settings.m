function [Settings,SubSetOutString] = instrument_settings(Instrument,Settings,SubSet)

%specific handling options for each inout observational dataset
%Corwin Wright, c.wright@bath.ac.uk, 01/March/2019


CoreVars = sampling_core_variables;

%identify and load information about the instrument
SubSetOutString = ''; %by default, no subsets.
switch Instrument
  case 'hirdls'
    Settings.ObsProperties.Path       = [CoreVars.MasterPath,'/tracks/HIRDLS/'];
    Settings.ObsProperties.FileString = 'hirdls';
    Settings.FineGrid                 = CoreVars.Hirdls.FineGrid;
  case 'saber'
    Settings.ObsProperties.Path       = [CoreVars.MasterPath,'/tracks/SABER/'];
    Settings.ObsProperties.FileString = 'saber';
    Settings.FineGrid                 = CoreVars.Saber.FineGrid;
  case 'cosmic'
    Settings.ObsProperties.Path       = [CoreVars.MasterPath,'/tracks/COSMIC/'];
    Settings.ObsProperties.FileString = 'cosmic';
    Settings.FineGrid                 = CoreVars.Cosmic.FineGrid;
  case 'airs'
    Settings.ObsProperties.Path       = [CoreVars.MasterPath,'/tracks/AIRS/'];
    Settings.ObsProperties.FileString = 'airs';
    Settings.FineGrid                 = CoreVars.Airs.FineGrid;
  case 'airs_qbo'
    Settings.ObsProperties.Path       = [CoreVars.MasterPath,'/tracks/AIRS_qbo/'];
    Settings.ObsProperties.FileString = 'qbo_airs';
    Settings.FineGrid                 = CoreVars.Airs.FineGrid;
  case 'qboz'
    Settings.ObsProperties.Path       = [CoreVars.MasterPath,'/tracks/QBOZ/'];
    Settings.ObsProperties.FileString = 'qbo_airs';
    Settings.FineGrid                 = CoreVars.Airs.FineGrid;
  case 'airs3d'
    Settings.ObsProperties.Path       = [CoreVars.MasterPath,'/tracks/AIRS_3D/'];
    Settings.ObsProperties.FileString = 'airs3d';
    Settings.FineGrid                 = CoreVars.Airs3D.FineGrid;
    SubSetOutString                   = ['_g',sprintf('%03d',SubSet)];
  case 'gisinger'
    Settings.ObsProperties.Path       = [CoreVars.MasterPath,'/tracks/gisinger/'];
    Settings.ObsProperties.FileString = 'airs3d';
    Settings.FineGrid                 = CoreVars.Airs3D.FineGrid;
    SubSetOutString                   = ['_g',sprintf('%03d',SubSet)];
  case 'airs_1d'
    Settings.ObsProperties.Path       = [CoreVars.MasterPath,'/tracks/AIRS_1D/'];
    Settings.ObsProperties.FileString = 'airs1d';
    Settings.FineGrid                 = CoreVars.Airs.FineGrid;
  case 'cosmic_petr'
    Settings.ObsProperties.Path       = [CoreVars.MasterPath,'/tracks/COSMICpetr/'];
    Settings.ObsProperties.FileString = 'cosmic_petr';
    Settings.FineGrid                 = CoreVars.Cosmic.FineGrid;    
  otherwise
    disp('Instrument not on specified list. Stopping')
    stop
end






end

