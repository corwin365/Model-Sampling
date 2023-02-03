function Weights = gong_channel_weights(Pressures)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%compute appropriately-weighted 3D AIRS brightness temperatures
%Corwin Wright, corwin.wright@trinity.oxon.org
%28/OCT/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%specify Gong et al channels, which we use
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Gong.Pressures = [2,2.5,3,4,7,10,20,30,40,60,80,100];
Gong.Levels    = zeros(numel(Gong.Pressures),14);
Gong.Levels(1,1)     = 74; %2.0 hPa
Gong.Levels(2,1)     = 75; %2.5 hPa
Gong.Levels(3,1)     = 76; %3.0 hPa
Gong.Levels(4,1)     = 77; %4.0 hPa
Gong.Levels(5,1)     = 78; %7.0 hPa
Gong.Levels(6,1)     = 79; %10.0 hPa
Gong.Levels(7,1:2)   = [81,82]; %20hPa
Gong.Levels(8,1:6)   = [102, 108, 114, 120, 125, 126 ] ;%30 hPa
Gong.Levels(9,1:7)   = [64, 88, 90, 94, 100, 106, 118]; %40 hPa
Gong.Levels(10,1:9)  = [66, 68, 70, 86, 87, 91, 93, 97, 130 ]; %60hPa
Gong.Levels(11,1:14) = [92, 98, 104, 105, 110, 111, 116, 117, 122, 123, 128, 129, 134, 140]; %80hPa
Gong.Levels(12,1:6)  = [132, 133, 138, 139, 149, 152 ]; %100 hPa

Gong.MinDetectable = sqrt(1e-3.*[3.78,3.72,3.63,3.66,3.88,4.62,2.14,0.98,0.83,0.66,0.50,0.67]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%now, for each pressure level we actually want, define a characteristic
%weighting function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AirsWeights = load([LocalDataDir,'/AIRS/airsweights.mat']);

Weighting = NaN(numel(Pressures),numel(AirsWeights.Pressures));


for iPrs=1:1:numel(Pressures)
  
  [~,levidx] = min(abs(Gong.Pressures-Pressures(iPrs)));
  Channels = Gong.Levels(levidx,:);
  Channels = Channels(Channels ~= 0);
  
  %take the mean of all WFs contributing
  Weighting(iPrs,:) = nanmean(AirsWeights.Weights(Channels,:),1);

end

Weights.Weights = Weighting;
Weights.PrsScale = AirsWeights.Pressures;

return

  