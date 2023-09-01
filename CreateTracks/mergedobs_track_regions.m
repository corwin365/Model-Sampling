clearvars
addpath(genpath('../'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%produce OIF data for COSMIC, HIRDLS, SABER and MLS
%
%divide into smaller geographic regions and temporal chunks, in order to 
%better fit very large models into memory
%
%
%Corwin Wright, c.wright@bath.ac.uk, 2023/08/26
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%dataset identifier
Settings.Instrument = 'limb_regions';

%geolocation - which data should we include?
%for all except HeightRange, we include any wholegranule including these
%for HeightRange, we will trim the granules in height to just this range
Settings.LatRange    = [-90,90];
Settings.LonRange    = [-180,180];
Settings.TimeRange   = datenum(2020,1,45);%:1:datenum(2020,1,20);
Settings.HeightScale = 20:0.1:80; %km

Settings.LatStep = 45;
Settings.LonStep = 45;
Settings.TimeStep = 3./24;

%path handling internal to routine
[~,CoreSettings] = sampling_core_v2(' ',' ',0,'GetSettings',true);
Settings.OutDir  = [CoreSettings.MasterPath,'/tracks/',Settings.Instrument,'/'];
clear CoreSettings


%individual instrument settings
%the instrument weighting functions are riddled with special cases, so we will handle them in-loop below, ugh. 
%Search for 'instrument weighting function shape' to edit these.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%GNSS
Settings.Inst.GNSS.ExtraVars = {'Azim'};

% %HIRDLS
% Settings.Inst.HIRDLS.ExtraVars = {};

%MLS
Settings.Inst.MLS.ExtraVars = {'LineOfSightAngle'};

%SABER
Settings.Inst.SABER.ExtraVars = {};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%extract lat,lon,z of each point, to allow model sampling
%store in daily files of a common format
%
%unlike the others, also store T, since it's downsampled.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for iDay=min(Settings.TimeRange):1:max(Settings.TimeRange);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %generate file name
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  disp('==========================')
  disp(datestr(iDay))

  %create a directory to put the data in, if it doesn't exist
  if exist(Settings.OutDir,'dir') ~= 7; mkdir(Settings.OutDir); end

  % % % %and name the output file
  % % % DayFile = [Settings.OutDir,'track_',Settings.Instrument,'_',num2str(iDay),'.mat'];
  % % % if exist(DayFile); continue; end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% import data
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  Instruments = fieldnames(Settings.Inst);
  clear Store
  for iInst=1:1:numel(Instruments)

    disp(['Loading ',Instruments{iInst}])
    % try; 
      Data =  get_limbsounders(iDay,Instruments{iInst},'HeightScale',Settings.HeightScale, ...
                                                       'AdditionalVars',Settings.Inst.(Instruments{iInst}).ExtraVars);
    % catch; continue; end
    if numel(Data.Temp) == 0; continue; end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %insturment orientation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if     strcmp(Instruments{iInst},'GNSS') %GNSS: do nothing - it's specified in the files and called 'Azim'
    elseif strcmp(Instruments{iInst},'MLS'); Data.Azim = Data.LineOfSightAngle;  %MLS: it's specified in the files and called 'LineOfSightAngle'
    elseif strcmp(Instruments{iInst},'HIRDLS') | strcmp(Instruments{iInst},'SABER')
      %we need to do some geometry to work out the LOS

      %then find the travel direction at each point
      LatScale = nanmean(Data.Lat,2);
      LonScale = nanmean(Data.Lat,2);
      Azim = azimuth(LatScale(1:end-1),LonScale(1:end-1),LatScale(2:end),LonScale(2:end))';
      Azim = [NaN,Azim]; %add back missing point due to the above shifting around
      Azim(1:2) = Azim(3); %should be very close

      switch Instruments{iInst}
        case 'HIRDLS'; Data.Azim = repmat(wrapTo180(360-((360-Azim)+180+47))',1,size(Data.Alt,2)); %HIRDLS
        case 'SABER';  Data.Azim = repmat(wrapTo180(360-((360-Azim)-90))',    1,size(Data.Alt,2)); %SABER
      end

      clear LatScale LonScale Azim
    end

    %vertical angle is always zero, because these are limb sounders (or it doesn't matter, for GNSS)
    Data.Vert = zeros(size(Data.Azim));

    %specify the instrument it came from
    Data.Inst = ones(size(Data.Azim)).*iInst;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %instrument weighting function shape
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    switch Instruments{iInst}
      case 'GNSS'

        ShapeX = 270; %270km along-LOS
        ShapeY = 1.5; %1.5km across-LOS
        ShapeZ = 1.5; %1.5km vertical width
        
        Data.X = single(ones(size(Data.Lat)).*ShapeX./4); %/2 for stdevs, additional /2 for width rather than distance from centre
        Data.Y = single(ones(size(Data.Lat)).*ShapeY./4);
        Data.Z = single(ones(size(Data.Lat)).*ShapeZ./4);

        clear ShapeX ShapeY ShapeZ

      case 'HIRDLS'

        ShapeX  = 200; %200km along-LOS
        ShapeY  = 20;  %20km across-LOS
        ShapeZ1 = 1;   %1km vertical width in s'sphere
        ShapeZ2 = 2;   %2km vertical width in m'sphere


        %the values are therefore
        Data.X = single(ones(size(Data.Lat)).*ShapeX./4); %/2 for stdevs, additional /2 for width rather than distance from centre
        Data.Y = single(ones(size(Data.Lat)).*ShapeY./4);

        %separate z in m'sphere and s'sphere
        Data.Z = single(ones(size(Data.Lat)));
        LowZ  = find(p2h(Data.Pres(:)) <= 60); Data.Z( LowZ) = ShapeZ1./4;
        HighZ = find(p2h(Data.Pres(:))  > 60); Data.Z(HighZ) = ShapeZ2./4;
        
        clear ShapeX ShapeY ShapeZ LowZ HighZ ShapeZ1 ShapeZ2

      case 'MLS'
        %approximate shape of blob
        ShapeX  = 165; %165km along-LOS
        ShapeY1  = 6;  %6km across-LOS below m'sphere
        ShapeY2  = 12;  %12km across-LOS below m'sphere

        %%% Calculate weightings based on resolutions:
        % X
        Data.X = single(ones(size(Data.Lat)).*ShapeX./4); %/2 for stdevs, additional /2 for width rather than distance from centre
        % Y
        %separate Y in t'sphere and above t'sphere (defined by pressure level 0.01hPa)
        Data.Y = single(ones(size(Data.Lat)));
        LowY  = find(Data.Pres  > 0.01); Data.Y( LowY) = ShapeY1./4;
        HighY = find(Data.Pres <= 0.01); Data.Y(HighY) = ShapeY2./4;

        % Z
        % Z is defined by table 3.22.1 in MLS v5-0_data_quality_document
        % which gives the following pressure values:
        given_pressures = [0.00022, 0.00046, 0.001, 0.01, 0.1, 0.316, 1, 3.6, 10,14.7, 31.6, 56.2, 100, 215, 261];
        % and their corresponding resolutions:
        given_horz_resolution = [12, 13,12 11, 6.4, 8.1, 6.8, 5.5, 4.1, 3.9, 3.6, 3.7, 4.6,3.5, 3.8];
    
        % convert to height
        given_heights = p2h(given_pressures);
        % interpolate to nearest neighbor for our pressures.
        ShapeZ = interp1(given_heights, given_horz_resolution, p2h(Data.Pres), 'nearest');
        Data.Z = ShapeZ./4;
        clear ShapeX ShapeY1 ShapeY2 LowY HighY given_horz_resolution given_pressures given_heights ShapeZ

      case 'SABER'

        ShapeX = 300;
        ShapeY = 50; 
        ShapeZ = 2;
        
        Data.X = single(ones(size(Data.Lat)).*ShapeX./4); %/2 for stdevs, additional /2 for width rather than distance from centre
        Data.Y = single(ones(size(Data.Lat)).*ShapeY./4);
        Data.Z = single(ones(size(Data.Lat)).*ShapeZ./4);

        clear ShapeX ShapeY ShapeZ
    end

    %store
    if ~exist('Store','var'); Store = Data; 
    else                      Store =  cat_struct(Store,Data,1);
    end

    clear Data
  end; clear iInst


  %clean up bad data

  Bad = find(Store.Temp <    0 | Store.Temp > 400 ...
           | Store.Lat  <  -90 | Store.Lat  > 90  ...
           | Store.Lon  < -180 | Store.Lon  > 180 ...
           | Store.Pres < 1e-5 | Store.Pres > 1200);
  Fields = fieldnames(Store);
  for iField=1:1:numel(Fields)
    F = Store.(Fields{iField});
    F(Bad) = NaN;
    Store.(Fields{iField}) = F;
  end; clear iField F Bad

stop

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %prepare the necessary information to reconstruct the granules from 1D data
  %and keep track of which instrument is which!
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

  x = 1:1:size(Store.Lat,1);  
  z = 1:1:size(Store.Lat,2);
  [z,x] = meshgrid(z,x);
  Recon.x = single(x(:));
  Recon.z = single(z(:));
  Recon.I = Store.Inst(:);
  Recon.Insts = Instruments;
  clear x z Instruments


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %prep data for output
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   


  %we don't want profiles to be split across regions. create a geolocation array that is just a single lat and lon for each profile
  RegionLat  = repmat(nanmean( Store.Lat,2),1,size( Store.Lat,2));
  RegionLon  = repmat(nanmean( Store.Lon,2),1,size( Store.Lon,2));
  RegionTime = repmat(nanmean(Store.Time,2),1,size(Store.Time,2));


  %put data into a line, reduce to singles to save space, and put into the appropriate format
  Track.Time       = Store.Time(:); %needs to be double  
  Track.Prs        = single(Store.Pres(:));
  Track.Lat        = single(Store.Lat(:));
  Track.Lon        = single(Store.Lon(:));  
  Track.T          = single(Store.Temp(:));
  Track.ViewAngleH = single(Store.Azim(:));
  Track.ViewAngleZ = single(Store.Vert(:));

  Weight.X = single(Store.X(:));
  Weight.Y = single(Store.Y(:));
  Weight.Z = single(Store.Z(:));

  RegionLat        = RegionLat(:);
  RegionLon        = RegionLon(:);
  RegionTime       = RegionTime(:);

  clear Store



  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  %now divide up into regions, and save
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

  MasterTrack  = Track;  clear Track
  MasterRecon  = Recon;  clear Recon
  MasterWeight = Weight; clear Weight


  SubSetCount = 0;
  for LonBand=min(Settings.LonRange):Settings.LonStep:max(Settings.LonRange)
    for LatBand=min(Settings.LatRange):Settings.LatStep:max(Settings.LatRange)
      for TimeBand=iDay:Settings.TimeStep:iDay+1

        InRange.Lon  = inrange(RegionLon,  LonBand+[0,Settings.LonStep ]);
        InRange.Lat  = inrange(RegionLat,  LatBand+[0,Settings.LatStep ]);
        InRange.Time = inrange(RegionTime,TimeBand+[0,Settings.TimeStep]);

        idx = intersect(intersect(InRange.Lon,InRange.Lat),InRange.Time);
        if numel(idx) == 0; continue; end

        %standard variables
        SubSetCount = SubSetCount+1;
        Track  = reduce_struct(MasterTrack, idx);
        Weight = reduce_struct(MasterWeight,idx);

        %reconstruction needs some manipulation to make the profiles indices increase sequentially
        Recon  = reduce_struct(MasterRecon, idx,{'Insts'});
        [~,~,Recon.x] = unique(Recon.x);

        %save
        OutFile = [Settings.OutDir,'track_',Settings.Instrument,'_',num2str(iDay),'_r',sprintf('%03d',SubSetCount),'.mat'];
        save(OutFile,'Track','Recon','Weight');
      end
    end
  end
  clear LonBand LatBand TimeBand


  %tidy up, then done
  clear Track DayFile

  disp([datestr(iDay),' complete'])
end; clear iDay



