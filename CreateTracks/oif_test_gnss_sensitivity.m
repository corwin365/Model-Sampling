clearvars
addpath(genpath('../'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%produce OIF data for COSMIC, HIRDLS, SABER and MLS
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
Settings.Instrument = 'gnss_htest';
[~,CoreSettings] = sampling_core_v3(' ',' ',0,'GetSettings',true);
Settings.OutDir  = [CoreSettings.MasterPath,'/tracks/',Settings.Instrument,'/'];
clear CoreSettings

%instrument settings
Settings.Instruments = {'GNSS'};

%geolocation - which data should we include?
Settings.HeightScale = 18:0.1:45; %height grid to sample on.
Settings.LatRange    = [-90,90];
Settings.LonRange    = [-180,180];
Settings.HeightScale = 15:0.5:50; %km

%dates to load geolocation from, and dates to sample from. These must line up
%precisely if both exist. If only one exists, the same dates will be used for both
Settings.Dates.Geolocation = datenum(2020,1,23);
Settings.Dates.Sampling    = datenum(2020,1,23);

%size of geographical subset regions. To produce daily global files, set to values > [360,180,24]
Settings.Subsets = [30,30,3]; %deglon, deglat, hours



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% instrument parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%GNSS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Params.GNSS.ViewAngle     = 'Azim'; 
Params.GNSS.WeightType    = '1dgauss';

%diverges here from the original OIF-limb routine in order to allow sensitivity testing
Test.GNSS.X = [0,10,-10,20,-20,50,-50,100,-100]+270;
Test.GNSS.Y = [0,0.1,0.2,0.5,-0.5,-0.2,-0.1]+1.5;
Test.GNSS.Z = [0,0.1,0.2,0.5,-0.5,-0.2,-0.1]+1.5;
Params.GNSS.WeightDetails = [NaN,NaN,NaN]; %will be replaced in-loop with values selected from the above

[Test.GNSS.X,Test.GNSS.Y,Test.GNSS.Z] = ndgrid(Test.GNSS.X,Test.GNSS.Y,Test.GNSS.Z);

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

    %load the observational data for this dat
    Data = get_limbsounders(Settings.Dates.Geolocation(iDay),   ...
                            Settings.Instruments{iInst} ,       ...
                            'HeightScale',Settings.HeightScale, ...
                            'FileSource', true,                 ...
                            'TimeHandling',3,                   ...
                            'AdditionalVars',ExtraVars,         ...
                            'DateWarning',false);

    if numel(Data.Alt) == 0;
      disp(['--> No data found for ',Settings.Instruments{iInst}]);
      clear Data ExtraVars
      continue; 
    end
    Data = rmfield(Data,{'Temp'}); %where we're going we don't need T
    Data.Prs = Data.Pres; Data = rmfield(Data,'Pres'); %rename variable for backwards compatibility

    % %finalise the viewing angle
    if ~isnumeric(Params.(Settings.Instruments{iInst}).ViewAngle); 
      %copy over the observational azimuths to the field 'ViewAngleH'
      %we're doing it this slightly awkward way in case the variable
      %is itself called ViewAngleH, which none currently are
      a = Data.(Params.(Settings.Instruments{iInst}).ViewAngle);
      Data = rmfield(Data,ExtraVars);
      Data.ViewAngleH = a;
    else
      %find the travel direction at each point
      LatScale = nanmean(Data.Lat,2);
      LonScale = nanmean(Data.Lat,2);
      Azim = azimuth(LatScale(1:end-1),LonScale(1:end-1),LatScale(2:end),LonScale(2:end))';
      Azim = [NaN,Azim]; %add back missing point due to the above shifting around
      Azim(1:2) = Azim(3); %should be very close to true for any realistic orbital pattern
      Data.ViewAngleH = repmat(wrapTo180(360-((360-Azim)+Params.(Settings.Instruments{iInst}).ViewAngle))',1,size(Data.Alt,2));
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


%DIVERGES HERE TO ADD THE MULTIPLE TESTING STATES   
for iTest=1:1:numel(Test.GNSS.X)
  Params.(Settings.Instruments{iInst}).WeightDetails(1) = Test.GNSS.X(iTest);
  Params.(Settings.Instruments{iInst}).WeightDetails(2) = Test.GNSS.Y(iTest);
  Params.(Settings.Instruments{iInst}).WeightDetails(3) = Test.GNSS.Z(iTest);

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

    if exist('Store','var'); xbase = max(Store.Data.x,[],'all'); else; xbase = 0; end

    x = (1:1:size(Data.Lat,1))+xbase;
    z = (1:1:size(Data.Lat,2));
    [Data.z,Data.x] = meshgrid(z,x);

    clear xbase

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
      Store.Data.OriginalFiles = cat(1,Store.Data.OriginalFiles,Data.OriginalFiles);

      %weight
      Store.Weight        = cat_struct(Store.Weight,Weight,1,{'Format'});
      Store.Weight.Format = cat(1,Store.Weight.Format,Weight.Format);


    end
    clear Recon Weight
    
end; clear iTest
clear Data 

%%%%%%%%%%%%%%%%%%%%%%%%%%%GOES BACK TO NORMAL FROM HERE

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
  RegionLat  = repmat(nanmedian(Store.Data.Lat,2),1,size( Store.Data.Lat, 2));
  RegionLon  = repmat(nanmedian(Store.Data.Lon,2),1,size( Store.Data.Lon, 2));
  RegionTime = repmat(nanmedian(Store.Data.Time,2),1,size(Store.Data.Time,2)) + Settings.Dates.Geolocation(iDay) - Settings.Dates.Sampling(iDay);

  %divide data up into geographic regions
  disp('--> Generating output and dividing into regions')
  SubSetCount = 0;
  for LonBand=min(Settings.LonRange):Settings.Subsets(1):max(Settings.LonRange)
    for LatBand=min(Settings.LatRange):Settings.Subsets(2):max(Settings.LatRange)
      for TimeBand=(Settings.Dates.Geolocation(iDay):Settings.Subsets(3)/24:Settings.Dates.Geolocation(iDay)+1) 

        %increment file counter
        SubSetCount = SubSetCount+1;        

        %find data in the box
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        InRange.Lon  = inrange(RegionLon,  LonBand+[0,Settings.Subsets(1)    ],2);
        InRange.Lat  = inrange(RegionLat,  LatBand+[0,Settings.Subsets(2)    ],2);
        InRange.Time = inrange(RegionTime,TimeBand+[0,Settings.Subsets(3)./24],2);
        idx = intersect(intersect(InRange.Lon,InRange.Lat),InRange.Time);
        if numel(idx) == 0; continue; end        

        %pull out the necessary data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Track  = reduce_struct(Store.Data,  idx,{'OriginalFiles'},0);
        Weight = reduce_struct(Store.Weight,idx,{},0);

        %split into:
        % Recon: ONLY profile reconstruction data, to put the profiles back together
        % Track: track data, for computing how to sample, and also metadata for postprocessing
        %we'll also save 'Params', as it's used in some weight calculations
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        Recon = struct();
        Fields = {'x','z'};
        for iF=1:1:numel(Fields)
          Recon.(Fields{iF}) = Track.(Fields{iF});
          Track = rmfield(Track,Fields{iF});
        end
        %recon needs to be rebased for each output file. This is as fine as we 
        %retain the original profile numbers in 'Meta' anyway
        [~,~,Recon.x] =unique(Recon.x); [~,~,Recon.z] =unique(Recon.z);

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
  disp(['--> ',datestr(Settings.Dates.Geolocation(iDay)),' complete; ',num2str(SubSetCount),' track files written'])
  clear Store SubSetCount


end; clear iDay
