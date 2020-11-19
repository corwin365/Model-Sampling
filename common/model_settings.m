function [Model,OldData,Settings] = model_settings(DayNumber,ModelType,Settings,ObsGrid)

%specific handling options for each model dataset
%Corwin Wright, c.wright@bath.ac.uk, 01/March/2019


switch ModelType  
  case 'cesm_ck'; 
    Model = load_cesm_ck(       ObsGrid);
  case 'cfsr';    
    Model = load_cfsr(          DayNumber);  
  case 'ec_fc';   
    Model = load_ecmwf_forecast(DayNumber,Settings.HoursAhead,Settings.MaxPrs);
  case 'ecmwf_issi';   
    Model = load_ecmwf_issi(ObsGrid);     
  case 'era5';    
    Model = load_era5(          DayNumber,Settings.MaxPrs); 
  case 'erai';    
    Model = load_erai(          DayNumber);     
  case 'jra55';   
    Model = load_jra55(         DayNumber);
  case 'jra55c';  
    Model = load_jra55c(        DayNumber);  
  case 'merra2';  
    Model = load_merra2(        DayNumber);    
  case 'um_fc';   
    Model = load_um_forecast(   DayNumber,ObsGrid);  
  case 'um_issi';   
    Model = load_um_issi(   ObsGrid);    
  otherwise
    disp('Model not on specified list. Stopping')
    OldData.ModelID = ''; %we need to set this for the return
    Error = 2;
    return
end

OldData = Model;

end

