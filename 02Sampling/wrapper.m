

%  %  %settings
Models = {'ec_fc'};
Instruments = {'airs3d'};
TimeScale = [datenum(2018,1,224:1:230)];
ForecastHours = 240;%[0:3:72,84:12:240];
NoClobber = 1;
Granules = [1,240];

if ~exist('Granules'); Granules = 0; else Granules = Granules(1):1:Granules(2); end
if ~exist('ForecastHours'); ForecastHours = 0;  end

disp(datestr(clock))


  %put a dummy oldfile in place
  %this is done inside the loop to tidy up memory - hitting some trouble switching RAs
  clear OldData
  OldData.ModelID = ''; OldData.Time = [NaN,NaN];
  
for iForecast=1:1:numel(ForecastHours)
  for iDay=1:1:numel(TimeScale)
    for iInst = 1:1:numel(Instruments)
    for iModel=1:1:numel(Models)
      for iGranule=numel(Granules):-1:1


        %do the analysis
        disp('====================================')      
        tic
        [Error,OldData] = sampling_core_singlethread(Instruments{iInst},Models{iModel},TimeScale(iDay),OldData,NoClobber,[],Granules(iGranule),ForecastHours(iForecast));
        
        toc
        disp('====================================')
        end
      end
    end
  end
end
disp(datestr(clock))

%  %  addpath('../03Reconstruction')
%  %  reconstruct_data_loop
%  %  reconstruct_data_loop_2d
%  %  reconstruct_data_loop_3d
