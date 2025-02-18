clearvars
addpath(genpath('../'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%produce OIF data for proposed ALICE instrument
%
%can divide the data into smaller geographic regions and temporal chunks,
% in order to better fit very large models into memory.
%
%Corwin Wright, c.wright@bath.ac.uk, 2023/11/07
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%path handling
Settings.Instrument = 'alice';
[~,CoreSettings] = sampling_core_v3(' ',' ',0,'GetSettings',true);
Settings.OutDir  = [CoreSettings.MasterPath,'/tracks/',Settings.Instrument,'/'];
clear CoreSettings

%instrument settings
Settings.Instruments = {'ALICE'}; %list of instruments to include

%geolocation - which data should we include?
% Settings.HeightScale = [0.,0.67,1.33,2,2.67,3.33,4,4.67,5.33,6,6.67,7.33,8,8.67,9.33,10.,10.68,11.33,12,12.66,13.32,14,14.66,15.34,16,16.67,17.33,18,18.67,19.33,20.,20.67,21.33,22,22.67,23.33,24,24.67,25.33,26,26.68,27.33,28,28.66,29.32,30.,30.66,31.34,32,32.67,33.33,34,34.67,35.33,36,36.67,37.33,38,38.67,39.33,40.,40.67,41.33,42,42.68,43.33,44,44.66,45.32,46,46.66,47.34,48,48.67,49.33,50.,50.67,51.33,52,52.67,53.33,54,54.67,55.33,56,56.67,57.33,58,58.68,59.33,60.,60.66,61.32,62,62.66,63.34,64,64.67,65.33,66,66.67,67.33,68,68.67,69.33,70.,70.67,71.33,72,72.67,73.33,74,74.68,75.33,76,76.66,77.32,78,78.66,79.34,80];
% Settings.Prs         = [1000	908.500000000000	825.400000000000	750	681.300000000000	619	562	511	464.200000000000	421.700000000000	383.100000000000	348.100000000000	316.200000000000	287.300000000000	261	237.100000000000	215	195.700000000000	177.800000000000	161.600000000000	147	133.400000000000	121.200000000000	110	100	90.8500000000000	82.5400000000000	75	68.1300000000000	61.9000000000000	56.2000000000000	51.1000000000000	46.4200000000000	42.1700000000000	38.3100000000000	34.8100000000000	31.6200000000000	28.7300000000000	26.1000000000000	23.7100000000000	21.5000000000000	19.5700000000000	17.7800000000000	16.1600000000000	14.7000000000000	13.3400000000000	12.1200000000000	11	10	9.08500000000000	8.25400000000000	7.50000000000000	6.81300000000000	6.19000000000000	5.62000000000000	5.11000000000000	4.64200000000000	4.21700000000000	3.83100000000000	3.48100000000000	3.16200000000000	2.87300000000000	2.61000000000000	2.37100000000000	2.15000000000000	1.95700000000000	1.77800000000000	1.61600000000000	1.47000000000000	1.33400000000000	1.21200000000000	1.10000000000000	1	0.908500000000000	0.825400000000000	0.750000000000000	0.681300000000000	0.619000000000000	0.562000000000000	0.511000000000000	0.464200000000000	0.421700000000000	0.383100000000000	0.348100000000000	0.316200000000000	0.287300000000000	0.261000000000000	0.237100000000000	0.215000000000000	0.195700000000000	0.177800000000000	0.161600000000000	0.147000000000000	0.133400000000000	0.121200000000000	0.110000000000000	0.100000000000000	0.0908500000000000	0.0825400000000000	0.0750000000000000	0.0681300000000000	0.0619000000000000	0.0562000000000000	0.0511000000000000	0.0464200000000000	0.0421700000000000	0.0383100000000000	0.0348100000000000	0.0316200000000000	0.0287300000000000	0.0261000000000000	0.0237100000000000	0.0215000000000000	0.0195700000000000	0.0177800000000000	0.0161600000000000	0.0147000000000000	0.0133400000000000	0.0121200000000000	0.0110000000000000	0.0100000000000000];
Settings.HeightScale = [8,8.67,9.33,10.,10.68,11.33,12,12.66,13.32,14,14.66,15.34,16,16.67,17.33,18,18.67,19.33,20.,20.67,21.33,22,22.67,23.33,24,24.67,25.33,26,26.68,27.33,28,28.66,29.32,30.,30.66,31.34,32,32.67,33.33,34,34.67,35.33,36,36.67,37.33,38,38.67,39.33,40.,40.67,41.33,42,42.68,43.33,44,44.66,45.32,46,46.66,47.34,48,48.67,49.33,50.,50.67,51.33,52,52.67,53.33,54,54.67,55.33,56,56.67,57.33,58,58.68,59.33,60.,60.66,61.32,62,62.66,63.34,64,64.67,65.33,66,66.67,67.33,68,68.67,69.33,70.,70.67,71.33,72,72.67,73.33,74,74.68,75.33,76,76.66,77.32,78,78.66,79.34,80];
Settings.Prs         = [316.200000000000	287.300000000000	261	237.100000000000	215	195.700000000000	177.800000000000	161.600000000000	147	133.400000000000	121.200000000000	110	100	90.8500000000000	82.5400000000000	75	68.1300000000000	61.9000000000000	56.2000000000000	51.1000000000000	46.4200000000000	42.1700000000000	38.3100000000000	34.8100000000000	31.6200000000000	28.7300000000000	26.1000000000000	23.7100000000000	21.5000000000000	19.5700000000000	17.7800000000000	16.1600000000000	14.7000000000000	13.3400000000000	12.1200000000000	11	10	9.08500000000000	8.25400000000000	7.50000000000000	6.81300000000000	6.19000000000000	5.62000000000000	5.11000000000000	4.64200000000000	4.21700000000000	3.83100000000000	3.48100000000000	3.16200000000000	2.87300000000000	2.61000000000000	2.37100000000000	2.15000000000000	1.95700000000000	1.77800000000000	1.61600000000000	1.47000000000000	1.33400000000000	1.21200000000000	1.10000000000000	1	0.908500000000000	0.825400000000000	0.750000000000000	0.681300000000000	0.619000000000000	0.562000000000000	0.511000000000000	0.464200000000000	0.421700000000000	0.383100000000000	0.348100000000000	0.316200000000000	0.287300000000000	0.261000000000000	0.237100000000000	0.215000000000000	0.195700000000000	0.177800000000000	0.161600000000000	0.147000000000000	0.133400000000000	0.121200000000000	0.110000000000000	0.100000000000000	0.0908500000000000	0.0825400000000000	0.0750000000000000	0.0681300000000000	0.0619000000000000	0.0562000000000000	0.0511000000000000	0.0464200000000000	0.0421700000000000	0.0383100000000000	0.0348100000000000	0.0316200000000000	0.0287300000000000	0.0261000000000000	0.0237100000000000	0.0215000000000000	0.0195700000000000	0.0177800000000000	0.0161600000000000	0.0147000000000000	0.0133400000000000	0.0121200000000000	0.0110000000000000	0.0100000000000000];

Settings.LatRange    = [-90,90];
Settings.LonRange    = [-180,180];

%dates to load geolocation from, and dates to sample from. These must line up
%precisely if both exist. If only one exists, the same dates will be used for both
Settings.Dates.Geolocation = datenum(2030,2,2); %not a typo
Settings.Dates.Sampling    = datenum(2020,1,20);

%size of geographical subset regions. To produce daily global files, set to values > [360,180,24]
Settings.Subsets = [60,30,3]; %deglon, deglat, hours

%when using 2d weighting functions, discard points where the value is less than this propertion of the max(abs())
Settings.MinVal = 0.01;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% instrument parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%each instrument has:
%  ViewAngle     - the viewing angle from north, either as:
%                  > a variable name to load, or 
%                  > a numerical offset from forward-along-track direction, in degrees
%  WeightType    - an instruction of how the sampling weights are specified, selected from:
%                 > '1dgauss'   - a 1D Gaussian of fixed size in the along-LOS, across-LOS and vertical direction
%                 > '1dgauss_z' - as above, but height-varying
%                 > '2d_field'  - file containing a 2D function for each point in x-z space
%  WeightDetails - varies for each weighttype, as follows:
%                 > '1dgauss'    - a three-element vector of the [aLOS,xLOS,z] FWHM in km
%                 > '1dgauss_z'  - a 2xm-element field named 'X' containing in each row a height and an aLOS FWHM [km], and
%                                - a 2xn-element field named 'Y' containing in each row a height and an xLOS FWHM [km], and
%                                - a 2xp-element field named 'Z' containing in each row a height and a vertical FWHM [km]
%                                - the function will interpolate between the specified values
%                 > '2d_field'   - 'File' - a path to the file containing the 2D weight matrices in the aLOS-z plane
%                                - 'Y'    - a 1D estimate of the FWHM in the xLOS direction




%ALICE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Params.ALICE.ViewAngle     = 180; 
Params.ALICE.WeightType    = '1dgauss';
Params.ALICE.WeightDetails = [100,10,1];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% check inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%sampling and geolocation time series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DatesExist = [0,0];
if isfield(Settings.Dates,'Geolocation'); DatesExist(1) = 1; end
if isfield(Settings.Dates,'Sampling'   ); DatesExist(2) = 1; end
if sum(DatesExist) == 0; error('No dates specified'); end
if sum(DatesExist) == 2;
  if numel(Settings.Dates.Geolocation) ~= numel(Settings.Dates.Sampling)
    error('Geolocation and sampling time series must have the same number of points')
  end
end
if sum(DatesExist) == 1;
  if DatesExist(1) == 1; Settings.Dates.Sampling = Settings.Dates.Geolocation;
  else;                  Settings.Dates.Geolocation = Settings.Dates.Sampling;
  end
end
clear DatesExist

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% main processing loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iDay=1:1:numel(Settings.Dates.Geolocation)

  if Settings.Dates.Geolocation(iDay) ~= Settings.Dates.Sampling(iDay); disp(['Producing OIF data for ', datestr(Settings.Dates.Geolocation(iDay)),' with ',datestr(Settings.Dates.Sampling(iDay)),' sampling dates'])
  else;                                                                 disp(['Producing OIF data for ', datestr(Settings.Dates.Geolocation(iDay))])
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% load and prep observational data
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  for iInst=1:1:numel(Settings.Instruments)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %load data, including horiz viewing angle
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %do we need to load the viewing angle data?
    if isnumeric(Params.(Settings.Instruments{iInst}).ViewAngle); 
      ExtraVars = {};
    else 
      ExtraVars = {Params.(Settings.Instruments{iInst}).ViewAngle};
    end

    %load the observational track file
    Grid = rCDF([LocalDataDir,'/corwin/alice/ALICE_800km_orbit_Jan_launch.nc']);

    %convert datetime to Matlab units
    y = floor(Grid.date./10000);
    m = floor(Grid.date./100)-y.*100;
    d = floor(Grid.date)-y.*10000 - m.*100;
    time = datenum(y,m,d)+Grid.time./24;
    clear y m d

    %select day of interest
    OnThisDay = find(floor(time) == Settings.Dates.Geolocation(iDay));
    Grid.Lat = Grid.Lat(:,OnThisDay); Grid.Lon = Grid.Lon(:,OnThisDay); time = time(OnThisDay);

    %put everything into a get_limbsounders-style struct
    Data.Lat  = repmat(Grid.Lat,1,1,numel(Settings.HeightScale));
    Data.Lon  = repmat(Grid.Lon,1,1,numel(Settings.HeightScale));
    Data.Time = repmat(time',32,1,numel(Settings.HeightScale));
    Data.Alt  = repmat(permute(Settings.HeightScale,[1,3,2]),32,size(Data.Lon,2),1);
    Data.Prs  = repmat(permute(Settings.Prs,[1,3,2]),32,size(Data.Lon,2),1);
    clear time Grid



    % %finalise the viewing angle
    if ~isnumeric(Params.(Settings.Instruments{iInst}).ViewAngle); 
      %copy over the observational azimuths to the field 'ViewAngleH'
      %we're doing it this slightly awkward way in case the variable
      %is itself called ViewAngleH, which none currently are
      a = Data.(Params.(Settings.Instruments{iInst}).ViewAngle);
      Data = rmfield(Data,ExtraVars);
      Data.ViewAngleH = a;
    else
      %find the travel direction at each point, using the centre track
      LatScale = squeeze(Data.Lat(16,:,1));
      LonScale = squeeze(Data.Lon(16,:,1));

      Azim = azimuth(LatScale(1:end-1),LonScale(1:end-1),LatScale(2:end),LonScale(2:end));
      Azim = [NaN,Azim]; %add back missing point due to the above shifting around
      Azim(1:2) = Azim(3); %should be very close to true for any realistic orbital pattern
      Data.ViewAngleH = repmat(wrapTo180(360-((360-Azim)+Params.(Settings.Instruments{iInst}).ViewAngle)),32,1,size(Data.Alt,3));
    end

    %add an identifier for the instrument, so we can merge them together
    Data.Inst = ones(size(Data.Lat)).*iInst;

    %replace geolocation dates with sampled dates
    if Settings.Dates.Geolocation(iDay) ~= Settings.Dates.Sampling(iDay);
      Data.Time = Data.Time - Settings.Dates.Geolocation(iDay) + Settings.Dates.Sampling(iDay);
    end    

    %vertical angle is always zero, because these are limb sounders (or it doesn't matter, for GNSS)
    Data.ViewAngleZ = zeros(size(Data.ViewAngleH));
    

    clear ExtraVars a 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %specify sensing geometry
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    switch Params.(Settings.Instruments{iInst}).WeightType;

      case '1dgauss'
        %same 1D gaussian set everywhere. Pass numerical values through to function.
        Weight.Format = repmat({'gaussian'},size(Data.Lat));
        Weight.X = single(ones(size(Data.Lat))).*Params.(Settings.Instruments{iInst}).WeightDetails(1);
        Weight.Y = single(ones(size(Data.Lat))).*Params.(Settings.Instruments{iInst}).WeightDetails(2);
        Weight.Z = single(ones(size(Data.Lat))).*Params.(Settings.Instruments{iInst}).WeightDetails(3);   

      case '1dgauss_z'
        %a height-varying set of 1D Gaussians. Pass numerical vales through to function.
        Weight.Format = repmat({'gaussian'},size(Data.Lat));
        for Dim={'X','Y','Z'};
          w = Params.(Settings.Instruments{iInst}).WeightDetails.(Dim{1});
          if size(w,1) == 1; Weight.(Dim{1}) = single(ones(size(Data.Lat))).*w(2);
          else;              Weight.(Dim{1}) = repmat(interp1(w(:,1),w(:,2),nanmean(Data.Alt,1)),size(Data.Lat,1),1);
          end
        end; clear Dim w

      case '2d_field'
        %2D (aLOS, z) field, with a fixed xLOS Gaussian FWHM. 
        %Make dummy fields for the values, and pass through the details in the 'Params' field.
        Weight.Format = repmat({['specified_2d_',Settings.Instruments{iInst}]},size(Data.Lat));

        %dummy values to keep weight arrays the right size when mixing instruments
        Weight.X = single(NaN(size(Data.Lat)));
        Weight.Y = single(NaN(size(Data.Lat)));
        Weight.Z = single(NaN(size(Data.Lat)));   
      otherwise
        error('Routine not currently configured for this geometry type')
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %clean up bad data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Bad = find(Data.Lat  <  -90 | Data.Lat  > 90  ...
             | Data.Lon  < -180 | Data.Lon  > 180 ...
             | Data.Prs  < 1e-5 | Data.Prs  > 1200);
    Fields = fieldnames(Data); 
    for iField=1:1:numel(Fields); 
      if strcmpi(Fields{iField},'OriginalFiles'); continue; end
      F = Data.(Fields{iField}); F(Bad) = NaN;  Data.(Fields{iField}) = F;
    end; 
    clear iField F Bad Fields

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %prepare information needed for reconstruction 
    % of the original profiles
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [Data.x,Data.y,Data.z] = ndgrid(1:1:size(Data.Lat,1),1:1:size(Data.Lat,2),1:1:size(Data.Lat,3)); %grid


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %glue this dataset to the others
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~exist('Store','var'); 
      Store.Data   = Data;  
      Store.Weight = Weight;
      Store.Recon.Insts = Settings.Instruments;
    else
      %data
      Data.SourceFile = Data.SourceFile + nanmax(Store.Data.SourceFile,[],'all'); %shift the SourceFile indices to account for the previous instruments
      Store.Data = cat_struct(Store.Data,Data,1,{'OriginalFiles'});
      Store.Data.OriginalFiles = cat(1,Store.Data.OriginalFiles,flatten(Data.OriginalFiles));

      %weight
      Store.Weight        = cat_struct(Store.Weight,Weight,1,{'Format'});
      Store.Weight.Format = cat(1,Store.Weight.Format,Weight.Format);


    end
    clear Data Recon Weight
    

    disp(['--> Imported observational grid for ',Settings.Instruments{iInst}])
    clear a ExtraVars LonScale LatScale Azim Data
  end; clear iInst 

  if ~exist('Store','var');
    disp('-->No data found for this day from any instrument, skipping to next day')
    continue
  end




  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% output to track files
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  %create a directory to put the data in, if it doesn't exist
  if exist(Settings.OutDir,'dir') ~= 7; mkdir(Settings.OutDir); end  

  %we don't want profiles to be split across regions. create a geolocation array that is just a single lat and lon for each profile
  RegionLat  = squeeze(nanmedian(Store.Data.Lat, [1,3]));
  RegionLon  = squeeze(nanmedian(Store.Data.Lon, [1,3]));
  RegionTime = squeeze(nanmedian(Store.Data.Time, [1,3]) + Settings.Dates.Geolocation(iDay) - Settings.Dates.Sampling(iDay));


  %divide data up into geographic regions
  disp('--> Generating output and dividing into regions')
  SubSetCount = 0;
  for LonBand=min(Settings.LonRange):Settings.Subsets(1):max(Settings.LonRange)
    for LatBand=min(Settings.LatRange):Settings.Subsets(2):max(Settings.LatRange)
      for TimeBand=(Settings.Dates.Geolocation(iDay):Settings.Subsets(3)/24:Settings.Dates.Geolocation(iDay)+1) 

        %increment file counter
        SubSetCount = SubSetCount+1      


        %find data in the box
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        InRange.Lon  = inrange(RegionLon,  LonBand+[0,Settings.Subsets(1)    ],2);
        InRange.Lat  = inrange(RegionLat,  LatBand+[0,Settings.Subsets(2)    ],2);
        InRange.Time = inrange(RegionTime,TimeBand+[0,Settings.Subsets(3)./24],2);
        idx = intersect(intersect(InRange.Lon,InRange.Lat),InRange.Time);
        if numel(idx) == 0; continue; end        

        %pull out the necessary data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Track  = reduce_struct(Store.Data,  idx,{'OriginalFiles'},2);
        Weight = reduce_struct(Store.Weight,idx,{},2);

        %split into:
        % Recon: ONLY profile reconstruction data, to put the profiles back together
        % Track: track data, for computing how to sample, and also metadata for postprocessing
        %we'll also save 'Params', as it's used in some weight calculations
        %
        %all arrays also need flattening down to 1D
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Recon = struct();
        [Recon.x,Recon.y,Recon.z] = meshgrid(1:1:size(Track.Lat,1),1:1:size(Track.Lat,2),1:1:size(Track.Lat,3));
        Recon.x = Recon.x(:);Recon.y = Recon.y(:); Recon.z = Recon.z(:);

        F = fieldnames(Track);
        for iF=1:1:numel(F)
          if     strcmp(F{iF},'OriginalFiles');         continue;
          elseif strcmp(F{iF},'x') | strcmp(F{iF},'y') | strcmp(F{iF},'z'); Track =rmfield(Track,F{iF});
          else                                                              Track.(F{iF}) = flatten(Track.(F{iF}));
          end
        end
        F = fieldnames(Weight);
        for iF=1:1:numel(F)
          Weight.(F{iF}) = flatten(Weight.(F{iF}));
        end        
      

        Track.PatternDay = Settings.Dates.Geolocation(iDay); %what day did the geolocation data come from?
        Track.InstNames = Settings.Instruments';


        %save track file
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        OutFile = [Settings.OutDir,'track_',Settings.Instrument,'_',num2str(Settings.Dates.Sampling(iDay)),'_r',sprintf('%03d',SubSetCount),'.mat'];
        save(OutFile,'Track','Recon','Weight','Params');%,'Meta');


      end; clear TimeBand
    end; clear LatBand
  end; clear LonBand
  clear InRange idx Track Weight Recon RegionLon RegionLat RegionTime OutFile 



  %tidy up for next loop
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  disp(['--> ',datestr(Settings.Dates.Geolocation(iDay)),' complete; up to ',num2str(SubSetCount),' track files written'])
  clear Store SubSetCount


end; clear iDay
