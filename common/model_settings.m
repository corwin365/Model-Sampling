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

%regular gridded models. These can use either gridded or
%scattered interpolants, so we don't need to check this
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
ModelInfo.dyamond_icon5km.Arguments = 'ObsGrid,WantedModel';
ModelInfo.dyamond_icon5km.WantedModel = 'icon5km';

%models that load as lists of points rather than as grids
%these can only be used with the ScatteredInt flag set
%upstream, so check this later
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%DYAMOND-WINTER 5km ICON run, scattered data
ModelInfo.dyamond_icon5km_scat.FuncName  = 'load_dyamond_common_scattered';
ModelInfo.dyamond_icon5km_scat.Arguments = 'ObsGrid,WantedModel';
WantedModel = 'icon5km';
ModelInfo.ScatteredFlag = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OK, call the model data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if strcmp(ModelType,'JUSTLISTING')

  %we just want a list of possible models - return this as the output
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Model = fieldnames(ModelInfo);
  return

else

  %otherwise, call up the data from the model we asked for
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %if we're using a scattered-output model, check we set the necessary flag
  if isfield(ModelInfo.(ModelType),'ScatteredFlag')
    if ModelInfo.(ModelType).ScatteredFlag == 1
      if Settings.ScatteredInt ~= 1;
        disp('Chosen model loads as a list of points rather than a grid')
        disp('''ScatteredInt'' flag must be set to true in main call to use this model')
        disp('Stopping')
        stop
      end
    end
  end

  %do the call
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %one day I intend to write this without the 'eval' statement (I started below but got stuck),
  %but right now I'm tight for time. One day...
  
  %select model
  TheModel = ModelInfo.(ModelType);

  %identify specific model for functions that handle multiple models
  if isfield(TheModel,'WantedModel'); WantedModel = TheModel.WantedModel; end  

  %generate the eval() statement
  TheCall = ['Model = ',TheModel.FuncName,'(',TheModel.Arguments,');'];
  
  %...and make the call
  eval(TheCall);


  %commented bit below is a partial first attempt at a non-eval()-based call
  % % % %   %create a function handle for
  % % % %   FunctionHandle = str2func(TheModel.FuncName);
  % % % %
  % % % %   %parse the list of arguments and parameterise (this is where I got stuck!)
  % % % %   Args = strsplit(TheModel.Arguments,',');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% done! tidy up and return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%copy over to OldData - this is used to avoid multiple duplicate loads.
OldData = Model;

return
