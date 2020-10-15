function [Model,OldData,Settings] = model_settings(DayNumber,ModelType,Settings,ObsGrid)

%specific handling options for each model dataset
%Corwin Wright, c.wright@bath.ac.uk, 01/March/2019


switch ModelType
  case 'cfsr';    Model = load_cfsr(          DayNumber);
  case 'jra55';   Model = load_jra55(         DayNumber);
  case 'jra55c';  Model = load_jra55c(        DayNumber);
  case 'erai';    Model = load_erai(          DayNumber);     
  case 'merra2';  Model = load_merra2(        DayNumber);    
  case 'era5';    Model = load_era5(          DayNumber,Settings.MaxPrs);   
  case 'ec_fc';   Model = load_ecmwf_forecast(DayNumber,Settings.HoursAhead,Settings.MaxPrs);
  case 'um_fc';   Model = load_um_forecast(   DayNumber,ObsGrid);    
  case 'cesm_ck'; Model = load_cesm_ck(       ObsGrid);
  otherwise
    disp('Model not on specified list. Stopping')
    OldData.ModelID = ''; %we need to set this for the return
    Error = 2;
    return
end

OldData = Model;

end

