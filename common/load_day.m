function [Data,Error] = load_day(DayNumber,Instrument,Model,OldData,TI)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%logic loop handlng how data is called
%
%Corwin Wright, c.wright@bath.ac.uk
%08/Jan/2018
%
%error conditions:
%0 - successful
%1 - unable to find sample file
%2 - instrument reading function not written
%3 - cannot find obs data for this day
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%get paths
CoreVars = sampling_core_variables;

if ~exist('TI'); TimeInterp = '_timeinterp';
else
  if TI == 0; TimeInterp = '';
  else        TimeInterp = '_timeinterp';
  end  
end


if strcmp(Model,'real') ~=1;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %sampled model datasets
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  FilePath = [CoreVars.MasterPath,'/reconstructed/',Instrument,'/',Model, ...
                                  '/sampled_',num2str(DayNumber),TimeInterp,'.mat'];
  if ~exist(FilePath); Data = []; Error = 1; return;  end
  Data = load(FilePath); Data = Data.Sampled_Data;
  Error = 0; return;
  
else
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %observational dataset - laod and format like the samples
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if ~exist('OldData'); OldData.Name = ' '; end
  
  if     strcmp(Instrument,'SABER') == 1;
    Data = extract_saber_data(DayNumber,CoreVars.Saber.Path,OldData,1);
    if Data.Error ~= 0; Data =[]; Error = 3; return;  end;
  elseif strcmp(Instrument,'COSMIC_full') == 1; 
    Data = extract_cosmic_data(DayNumber,CoreVars.Cosmic.Path,OldData,1);
    if Data.Error ~= 0; Data =[]; Error = 3; return;  end;
  elseif strcmp(Instrument,'COSMIC') == 1;
    Data = extract_cosmic_data_fast(DayNumber,1);
    if Data.Error ~= 0; Data =[]; Error = 3; return;  end;    
  elseif strcmp(Instrument,'HIRDLS') == 1;
    Data = extract_hirdls_data(DayNumber,CoreVars.Hirdls.Path);
    if Data.Error ~= 0; Data =[]; Error = 3; return;    end;
  elseif strcmp(Instrument,'AIRS') == 1;
    stop
    % % %
    % % %     %AIRS-L1
    % % %     if Data.Error ~= 0; Data =[]; Error = 3; return;    end;
    % % %     stop
  else
    Data = [];Error = 2; return
  end
  
  
  if exist('OldData'); Data.OldData = OldData; end
  
  %success!
  Error = 0; return
  
  
end

