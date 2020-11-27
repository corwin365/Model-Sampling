% % function Model = load_cesm_ck(ObsGrid)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load CESM model data from Chris, and put into common analysis format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%possible errors:
%0: success
%2. file not found


%get core variables - needed for model data path
CoreVars = sampling_core_variables;
% CoreVars.CESM_CK.Path   = [LocalDataDir,'/CESM/'];
CoreVars.CESM_CK.Path   = [LocalDataDir,'/corwin/issi/chris/'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load the data for this input grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%define timestep grid
Step = 1./24./60.*15;
First = min(ObsGrid.Track.Time(:));
Last  = max(ObsGrid.Track.Time(:));

%and pressure grid
MaxPrs = 10.^max(ObsGrid.Track.Prs(:)).*1.2;
MinPrs = 10.^min(ObsGrid.Track.Prs(:))./1.2;


%round 'first' and 'last' time to the nearest 15-min step, rounding away from the sample
First = floor(First./Step).*Step;
Last  = ceil(  Last./Step).*Step;

for Time=First:Step:Last;

  
  %identify file
  [y,M,d,h,m,~] = datevec(Time);
% %   FileName = [CoreVars.CESM_CK.Path, ...
% %               'wrfout_d01_',sprintf('%04d',y),'-',sprintf('%02d',M),'-',sprintf('%02d',d), ...
% %                         '_',sprintf('%02d',h),'-',sprintf('%02d',m),'-00.nc'];
                 
  FileName = [CoreVars.CESM_CK.Path, ...
              'wrfout_d01_',sprintf('%04d',y),'-',sprintf('%02d',M),'-',sprintf('%02d',d), ...
                        '_',sprintf('%02d',h),'-',sprintf('%02d',m),'-00.mat'];
%   FileName = 'C:\Data\corwin\issi\chris\wrfout_d01_2010-10-08_12-00-00.mat' %temporary override for local testing                           
    
  %load file
  if ~exist(FileName,'file');
    Model.Error = 2;
    return
  end
%   Data = cjw_readnetCDF(FileName,1);  
  Data = load(FileName);  
  
% % % %   %convert units, as instructed by CK
% % % %   P = Data.PB + Data.P; P = P .* 0.01;
% % % % 
% % % %   %make pressure easier to work with
% % % %   Prs = squeeze(nanmean(P,[1,2]));%average pressure over horizontal plane (bad near surface, fine in s'sphere)
% % % %   
% % % %   %drop unwanted pressure levels
% % % %   Good = find(Prs <= MaxPrs & Prs >= MinPrs);
% % % %   Data.T  = Data.T( :,:,Good);
% % % %   P       = P(:,:,Good);
% % % %   Prs     = Prs(Good);
% % % %   
% % % %   %offset described by Chris Kruse
% % % %   Data.T = Data.T + 300;
% % % %   
% % % %   %convert theta to T
% % % %   T = single(Data.T)./((1000./P).^0.2896);


  %pull out and reformat data
  if Time == First;
    AllData.T    = Data.T;
    AllData.Prs  = Data.Prs;
    AllData.Time = Time;
    AllData.Lat  = Data.Lat;
    AllData.Lon  = Data.Lon';
  else
    AllData.T    = cat(4,AllData.T,Data.T);
    AllData.Time = cat(1,AllData.Time,Time);
    
  end
  
  clear Data Good Prs T FileName y M d h m 

end; clear Time First Last Step

% % % % %we need 1d lat and lon. So, we need to reinterpolate the data to a regular lat/lon grid
% % % % %oversample in both directions, to be safe
% % % % MinLat = max([min(ObsGrid.Track.Lat(:)),min(AllData.Lat(:))]);
% % % % MaxLat = min([max(ObsGrid.Track.Lat(:)),max(AllData.Lat(:))]);
% % % % MinLon = max([min(ObsGrid.Track.Lon(:)),min(AllData.Lon(:))]);
% % % % MaxLon = min([max(ObsGrid.Track.Lon(:)),max(AllData.Lon(:))]);
% % % % 
% % % % Lat = MinLat:0.05:MaxLat;
% % % % Lon = MinLon:0.05:MaxLon;
% % % % 
% % % % %  %  %  Lat = min(ObsGrid.Track.Lat(:))-1 : 0.08: max(ObsGrid.Track.Lat(:))+1;
% % % % %  %  %  Lon = min(ObsGrid.Track.Lon(:))-1 : 0.08: max(ObsGrid.Track.Lon(:))+1;
% % % % [xi,yi] = meshgrid(Lon,Lat);
% % % % 
% % % % %create interpolant object
% % % % F = scatteredInterpolant(double(flatten(AllData.Lon)), ...
% % % %                          double(flatten(AllData.Lat)), ...
% % % %                          double(flatten(AllData.T(:,:,1,1))));
% % % % 
% % % % %create storage array
% % % % sz = size(AllData.T);
% % % % T2 = NaN([size(xi),sz(3),sz(4)]);
% % % % 
% % % % %regrid
% % % % for iTime=1:1:sz(4);
% % % %   for iLevel=1:1:sz(3);
% % % %     F.Values = double(flatten(AllData.T(:,:,iLevel,iTime)));
% % % %     T2(:,:,iLevel,iTime) = F(double(xi),double(yi));
% % % %   end
% % % % end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%reformat for analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%output format:
%struct called Model
%containing fields:
%Lon   - 1d. Runs from -180 to +180
%Lat   - 1d
%Time  - 1d
%Prs   - 1d
%T     - 4d, time x lon x lat x pressure


% stick stuff in a struct
Model.Lon  = AllData.Lon;
Model.Lat  = AllData.Lat;
Model.Time = AllData.Time;
Model.T    = double(permute(AllData.T,[4,2,1,3]));
Model.Prs  = AllData.Prs;



%success!
Model.Error = 0;
return


