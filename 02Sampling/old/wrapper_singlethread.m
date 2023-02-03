clear all

%  %settings
Models = {'um_fc'};
Instruments = {'airs3d'};
TimeScale = 737284;
Granules = [3,3];
NoClobber = 1;
ForecastHours = 0;

if ~exist('Granules'); Granules = 0; else Granules = Granules(1):1:Granules(2); end

disp(datestr(clock))


  %put a dummy oldfile in place
  %this is done inside the loop to tidy up memory - hitting some trouble switching RAs
  clear OldData
  OldData.ModelID = ''; OldData.Time = [NaN,NaN];
  

  for iDay=1:1:numel(TimeScale)
    for iInst = 1:1:numel(Instruments)
    for iModel=1:1:numel(Models)
      for iGranule=1:1:numel(Granules)

        %do the analysis
        disp('====================================')      
        tic
        [Error,OldData] = sampling_core_singlethread(Instruments{iInst},Models{iModel},TimeScale(iDay),OldData,NoClobber,[],Granules(iGranule),ForecastHours);
        
        toc
        disp('====================================')
      end
    end
  end
end
disp(datestr(clock))

%  %  addpath('../03Reconstruction')
%  %  reconstruct_data_loop
%  %  reconstruct_data_loop_2d
%  %  reconstruct_data_loop_3d
