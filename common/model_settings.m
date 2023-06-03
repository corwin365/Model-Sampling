function [Model,OldData,Settings] = model_settings(DayNumber,ModelType,Settings,ObsGrid)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%specific handling options for each model dataset
%Corwin Wright, c.wright@bath.ac.uk, 01/March/2019
%heavily rewritten 13/October/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%produce a struct containing details of the models, the functions they are called by, and the arguments those functions need
%'JUSTLISTING' will be an override argument for just getting a list of possible models, so don't use this as a model name!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create the base struct, then we can start filling it
ModelInfo = struct;

%CESM data produced by Chris Kruse for ISSI 2019
ModelInfo.cesm_ck.FuncName  = 'load_cesm_ck';
ModelInfo.cesm_ck.Arguments = 'ObsGrid';

%CFSR reanalysis
ModelInfo.cfsr.FuncName  = 'load_cfsr';
ModelInfo.cfsr.Arguments = 'DayNumber';

%ECMWF Forecasts
ModelInfo.ec_fc.FuncName  = 'load_ecmwf_forecast';
ModelInfo.ec_fc.Arguments = 'DayNumber,Settings.HoursAhead,Settings.MaxPrs';

%ECMWF runs produced by Inna for ISSI 2019
ModelInfo.ecmwf_issi.FuncName  = 'load_ecmwf_issi';
ModelInfo.ecmwf_issi.Arguments = 'ObsGrid';

%ERA5 reanalysis
ModelInfo.era5.FuncName  = 'load_era5';
ModelInfo.era5.Arguments = 'DayNumber,Settings.MaxPrs,Settings.MinPrs';

%ERA-Interim reanalysis
ModelInfo.erai.FuncName  = 'load_erai';
ModelInfo.erai.Arguments = 'DayNumber';

%1.4km IFS
ModelInfo.ifs_1km.FuncName  = 'load_1kmIFS';
ModelInfo.ifs_1km.Arguments = 'ObsGrid';

%1.4km IFS grid with a fake wave replacing the data
ModelInfo.ifs_1km_fakewave.FuncName  = 'load_1kmIFS_fakewave';
ModelInfo.ifs_1km_fakewave.Arguments = 'ObsGrid';

%JRA55 reanalysis
ModelInfo.jra55.FuncName  = 'load_jra55';
ModelInfo.jra55.Arguments = 'DayNumber';

%JRA55C reanalysis
ModelInfo.jra55c.FuncName  = 'load_jra55c';
ModelInfo.jra55c.Arguments = 'DayNumber';

%MERRA2 reanalysis
ModelInfo.merra2.FuncName  = 'load_merra2';
ModelInfo.merra2.Arguments = 'DayNumber';

%UM runs produced by Annelize for ISSI 2019
ModelInfo.um_issi.FuncName  = 'load_um_issi';
ModelInfo.um_issi.Arguments = 'ObsGrid';

%DYAMOND-WINTER 5km UM runs
ModelInfo.dyamond_um5k.FuncName  = 'load_dyamond_um5ktest';
ModelInfo.dyamond_um5k.Arguments = 'ObsGrid';

%DYAMOND-WINTER 5km ICON run
ModelInfo.dyamond_icon5km.FuncName  = 'load_dyamond_common';
ModelInfo.dyamond_icon5km.Arguments = 'ObsGrid,WantedModel,FixedPFlag';
ModelInfo.dyamond_icon5km.WantedModel = 'icon5km';
ModelInfo.dyamond_icon5km.FixedPFlag  = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OK, call the model data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if strcmp(ModelType,'JUSTLISTING')
  %we just want a list of possible models - return this as the output
  Model = fieldnames(ModelInfo);
  return
else
  %otherwise, call up the data from the model we asked for
  %one day I intend to write this without the 'eval' statement (I started below but got stuck),
  %but right now I'm tight for time. One day...

  %identify model
  TheModel = ModelInfo.(ModelType);  

  %pull out 'WantedModel' and 'FixedPFlag' for the right model if requested
  if isfield(TheModel,'WantedModel'); WantedModel = TheModel.WantedModel; end
  if isfield(TheModel, 'FixedPFlag'); FixedPFlag  = TheModel.FixedPFlag;  end

  %generate eval() statement...
  TheCall = ['Model = ',TheModel.FuncName,'(',TheModel.Arguments,');'];

  %run it
  eval(TheCall);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% done! tidy up and return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%copy over to OldData - this is used to avoid multiple duplicate loads.
OldData = Model;

return
