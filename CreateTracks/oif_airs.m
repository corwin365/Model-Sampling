clearvars
addpath(genpath('../'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%produce OIF data for 3D AIRS. No significant changes from previous version,
%but significantly cleaned up to support any future changes
%
%Corwin Wright, c.wright@bath.ac.uk, 2023/11/07
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%path handling
Settings.Instrument = 'AIRS3D';
[~,CoreSettings] = sampling_core_v3(' ',' ',0,'GetSettings',true);
Settings.OutDir  = [CoreSettings.MasterPath,'/tracks/',Settings.Instrument,'/'];
clear CoreSettings

%geolocation - which data should we include?
Settings.LatRange    = [-90,90];
Settings.LonRange    = [-180,180];
Settings.HeightRange = [20,60]; %km

%dates to load geolocation from, and dates to sample from. These must line up
%precisely if both exist. If only one exists, the same dates will be used for both
Settings.Dates.Geolocation = datenum(2020,1,1);%20:1:60);
Settings.Dates.Sampling    = datenum(2020,1,1);%20:1:60);

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

  disp('=======================================================')
  if Settings.Dates.Geolocation(iDay) ~= Settings.Dates.Sampling(iDay); disp(['Producing OIF data for ', datestr(Settings.Dates.Geolocation(iDay)),' with ',datestr(Settings.Dates.Sampling(iDay)),' sampling dates'])
  else;                                                                 disp(['Producing OIF data for ', datestr(Settings.Dates.Geolocation(iDay))])
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% load and prep observational data
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  k= 0; %just a count for display to the user at the end
  for iGranule=1:1:240;
    
    OutFile = [Settings.OutDir,'track_',Settings.Instrument,'_',num2str(Settings.Dates.Sampling(iDay)),'_g',sprintf('%03d',iGranule),'.mat'];
    if exist(OutFile,'file'); continue;end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %load data, including horiz viewing angle
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %load the observational data for this day
    [Airs,~,Error] = prep_airs_3d(Settings.Dates.Sampling(iDay),iGranule,'LoadOnly',true);
    if Error ~= 0; clear Airs Error; continue; end %failed to load data
    clear Error

    %convert filename from path to just name
    [~,a,c] = fileparts(Airs.Source);
    Airs.Source = [a,c];
    clear a c

    %extract the fields we need
    Fields = {'l1_time','l1_lat','l1_lon','ret_z','Source'};
    Data = struct();
    for iF=1:1:numel(Fields); Data.(Fields{iF}) = Airs.(Fields{iF}); end
    clear Airs Fields iF 

    %repmat the fields to all match
    sz = [size(Data.l1_lon),numel(Data.ret_z)];
    Data.l1_lon  = repmat(Data.l1_lon, [1,1,sz(3)]);
    Data.l1_lat  = repmat(Data.l1_lat, [1,1,sz(3)]);   
    Data.l1_time = repmat(Data.l1_time,[1,1,sz(3)]); 
    Data.ret_z   = repmat(permute(Data.ret_z,[2,3,1]),[sz(1:2),1]);

    %create tracking information back to the source data
    Data.SourceProf = reshape(1:1:prod(sz),sz);
    Data.SourceFile = ones(sz);
    

    %create an approximate pressure field
    Data.Prs = h2p(Data.ret_z);

    %the vertical viewing angle varies by distance off-axis geometrically
    %angle is defined as a/c/w from nadir
    ViewAngleZ = repmat(permute((1:1:sz(1)) - (90./2),[2,1]),[1,sz(2:3)]); %rows off-centre
    ViewAngleZ(1:45,:,:) = ViewAngleZ(1:45,:,:)-1;
    ViewAngleZ = ViewAngleZ ./ max(abs(ViewAngleZ)); %normalised
    ViewAngleZ = ViewAngleZ .* 49.5; %degrees
    Data.ViewAngleZ = ViewAngleZ; clear ViewAngleZ


    %the horizontal viewing angle depends on the satellite scan track location
    %we need to use a row-based approach, otherwise we'll get a jump at the end of each height level
    LatScale = Data.l1_lat; NextLat = circshift(LatScale,[0,-1,0]); LatScale(:,end,:) = NaN; NextLat( :,end,:) = NaN;
    LonScale = Data.l1_lon; NextLon = circshift(LonScale,[0,-1,0]); LonScale(:,end,:) = NaN; NextLon( :,end,:) = NaN;
     
    %then find the travel direction at each point
    %and replace the point we removed with the next value, which should be very
    %close to being correct
    Azimuth = azimuth(LatScale,LonScale,NextLat,NextLon);    
    Azimuth(:,end,:) = Azimuth(:,end-1,:);
    
    %finally, convert to a viewing angle
    ViewAngleH = wrapTo180(360-((360-Azimuth)))+45; %degrees c/w from N
    
    %smooth over minor jumps caused by finite precision on the geolocation, and store
    %the nature of the data means large discontinuities don't exist, so this is safe
    Data.ViewAngleH = smoothn(ViewAngleH,[5,5,1]);
    clear LatScale LonScale NextLat NextLon Azimuth ViewAngleH

    %finally, do some renaming
    Data.Lat  = Data.l1_lat;  Data = rmfield(Data,'l1_lat');
    Data.Lon  = Data.l1_lon;  Data = rmfield(Data,'l1_lon');
    Data.Time = Data.l1_time; Data = rmfield(Data,'l1_time');
    Data.Alt  = Data.ret_z;   Data = rmfield(Data,'ret_z');
    Data.OriginalFiles = {Data.Source}; Data = rmfield(Data,'Source');

    %and replace geolocation dates with sampled dates
    if Settings.Dates.Geolocation(iDay) ~= Settings.Dates.Sampling(iDay);
      Data.Time = Data.Time - Settings.Dates.Geolocation(iDay) + Settings.Dates.Sampling(iDay);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %clean up bad data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %sometimes the files are missing height scales. If so, stick in the standard one for these data.
    if sum(Data.Alt(:)) == 0; Data.Alt = repmat(permute([0;3;6;9;12;15;18;21;24;27;30;33;36;39;42;45;48;51;54;57;60;65;70;75;80;85;90],[2,3,1]),[sz(1:2),1]);end

    %sometimes the data contain incorrect jumps in time. The granule should all be within a six-minute range
    %by definition, so if the data is more than 15min away from the median then remove and interpolate
    Data.Time(find(abs(Data.Time-median(Data.Time,'all')).*60.*24 > 15)) = NaN; 
    for iLev=1:1:sz(3); Data.Time(:,:,iLev) = inpaint_nans(Data.Time(:,:,iLev)); end; clear iLev

    %discard unwanted heights
    Data = reduce_struct(Data,inrange(nanmean(Data.Alt,[1,2]),(Settings.HeightRange)),{'OriginalFiles'},3);
    clear InHeightRange

    %remove unphysical error values 
    Bad = find(Data.Lat  <  -90 | Data.Lat  > 90  ...
             | Data.Lon  < -180 | Data.Lon  > 180 ...
             | Data.Prs  < 1e-5 | Data.Prs  > 1200);
    Fields = fieldnames(Data);
    for iField=1:1:numel(Fields);
      if strcmpi(Fields{iField},'OriginalFiles'); continue; end
      F = Data.(Fields{iField}); F(Bad) = NaN;  Data.(Fields{iField}) = F;
    end;
    clear iField F Bad Fields sz

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %specify sensing geometry
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %tell the sampler what kind of data it's getting
    Weight.Format = repmat({'gaussian'},size(Data.Lat));

    %the  horizontal weighting functions are approximately...
    sz = size(Data.Lat);    
    Weight.X = ones(sz).*13.5./(2.*2.355); %along-LOS. the factor of 2.355converts from FWHM to stdev, the factor of 2 from stdev to the stdev/2 used internally
    Weight.Y = ones(sz).*13.5./(2.*2.355); %across-LOS

    %vertical comes from the geolocation and whether it's day or night
    %we'll compute a single day/night profile value for the granule, as this is simpler and it shouldn't vary much,
    %and then use either the day or night value for each point as appropriate
    [latmean,~] = meanm(Data.Lat(:),Data.Lon(:));
    z = squeeze(nanmean(Data.Alt,[1,2]));
    RDay   = airs_resolution(0,date2doy(Settings.Dates.Geolocation(iDay)),latmean,z,0)./(2.*2.355);
    RNight = airs_resolution(1,date2doy(Settings.Dates.Geolocation(iDay)),latmean,z,0)./(2.*2.355);

    IsDay = which_airs_retrieval(Data.Lon,Data.Lat,Data.Time,1);
    Weight.Z = Weight.X.*NaN;
    Weight.Z(IsDay == 1) = interp1(z,  RDay,Data.Alt(IsDay == 1)); 
    Weight.Z(IsDay == 0) = interp1(z,RNight,Data.Alt(IsDay == 0));
    clear latmean z IsDay RDay RNight sz

    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %prepare information needed for reconstruction
    % of the original profiles
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [Data.x,Data.y,Data.z] = ndgrid(1:1:size(Data.Lat,1),1:1:size(Data.Lat,2),1:1:size(Data.Lat,3)); %grid

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %store data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %this segment of code has the functionality to append multiple granules
    %this currently does nothing, as the surrounding code is only working on
    %a single granule at the time, but the additional runtime is less than 
    %random jitter on the total time taken, so I've left it in place to support
    %future revision possibilities

    if ~exist('Store','var');
      Store.Data   = Data;
      Store.Weight = Weight;
      Store.Recon.Insts = Settings.Instrument;
    else
      %data
      Data.SourceFile = Data.SourceFile + nanmax(Store.Data.SourceFile,[],'all'); %shift the SourceFile indices to account for the previous instruments
      Store.Data = cat_struct(Store.Data,Data,2,{'OriginalFiles'});
      Store.Data.OriginalFiles = cat(2,Store.Data.OriginalFiles,Data.OriginalFiles);

      %weight
      Store.Weight        = cat_struct(Store.Weight,Weight,1,{'Format'});
      Store.Weight.Format = cat(2,Store.Weight.Format,Weight.Format);
    end
    clear Data  Weight


    disp(['--> Imported observational grid for granule ',num2str(iGranule)])
    clear a ExtraVars LonScale LatScale Azim Data


  if ~exist('Store','var');
    disp('-->No data found for this day, skipping to next day')
    continue
  end




  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% output to track files
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %create a directory to put the data in, if it doesn't exist
  if exist(Settings.OutDir,'dir') ~= 7; mkdir(Settings.OutDir); end  

  disp('--> Generating output')

  %pull out the necessary data, and:
  %1. flatten it into a 1D vector
  %2. if it's not time, make it a single
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  Track  = Store.Data;
  Weight = Store.Weight;


  sz = size(Track.Lat);

  f = fieldnames(Track);
  for iF=1:1:numel(f)
    if isequal(sz,size(Track.(f{iF}))); 
      Track.(f{iF}) = flatten(Track.(f{iF}));
      if ~strcmpi(f{iF},'Time');  Track.(f{iF}) = single(Track.(f{iF})); end
    end
  end; clear iF f

  f = fieldnames(Weight);
  for iF=1:1:numel(f); 
    Weight.(f{iF}) = flatten(Weight.(f{iF})); 
    if strcmpi(class(Weight.(f{iF})),'double'); Weight.(f{iF}) = single(Weight.(f{iF})); end
  end; clear f iF


  %split into:
  % Recon: ONLY profile reconstruction data, to put the profiles back together
  % Track: track data, for computing how to sample, and also metadata for postprocessing
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  Recon = struct();
  Fields = {'x','y','z'};
  for iF=1:1:numel(Fields)
    f = Track.(Fields{iF});
    Recon.(Fields{iF}) = f(:); clear f
    Track = rmfield(Track,Fields{iF});
  end; clear iF Fields

  Track.PatternDay = Settings.Dates.Geolocation(iDay); %what day did the geolocation data come from?
  Track.InstNames = Settings.Instrument;



  %save track file
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    save(OutFile,'Track','Recon','Weight');
  disp('--> Track file written')
  k = k+1;

  clear Track Recon Weight sz OutFile Store 


  end; clear iGranule


  %tidy up for next loop
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  disp(['----> ',datestr(Settings.Dates.Geolocation(iDay)),' complete; ',num2str(k),' track files written <----'])


end; clear iDay
