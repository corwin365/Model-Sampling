
  function [Fine,W] = gaussian_blob(Sample,ObsGrid,Settings)


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %compute fine grid and weight for blobs defined by 1D Gaussians in
  %three dimensions
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %produce the fine sampling grid
  %this is the grid the actual sampling will be done on
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %declare a struct to keep the grids in
  Fine = struct();

  %creating the grid in the horizontal is easy
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  Fine.X = -Settings.BlobScale*Sample.WeightX : Settings.FineGrid(1) : Settings.BlobScale*Sample.WeightX;
  Fine.Y = -Settings.BlobScale*Sample.WeightY : Settings.FineGrid(2) : Settings.BlobScale*Sample.WeightY;

  %for the vertical, we're working in pressure co-ordinates, so this is fiddly
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %most instruments are a simple Gaussian in height
  %these are specified in the track generator functions as having *POSITIVE* Sample.WeightZ values
  %so create a height window accordingly
  if Sample.WeightZ >= 0
    %first, find approx height of the measurement-centre pressure level
    z = p2h(10.^Sample.Prs);
    %next, find height above and below this that we want, and convert to log-P
    z = (Settings.BlobScale*[-1,1].*Sample.WeightZ)+z;
    z = log10(h2p(z));
  else
    %some instruments are more complex than just a simple Gaussian
    %these are specified with a *NEGATIVE* Sample.WeightZ
    %the (integer) values of this tells us which of a pre-defined set of functions to use
    %from those stored in the source file

    %identify which channel we want.
    Channel = struct();
    Channel.ID = abs(Sample.WeightZ);

    Channel.Prs = log10(ObsGrid.Weight.ZFuncs.PrsScale);
    Channel.W   = ObsGrid.Weight.ZFuncs.Weights(Channel.ID,:); %keep for later

    %discard parts that contribute very little
    Channel.W(Channel.W < Settings.MinZContrib.*sum(Channel.W(:)),'omitnan') = 0;

    %set NaNs to zero
    Channel.W(isnan(Channel.W)) = 0;

    %remove low altitudes
    GoodPrs = find(Channel.Prs < log10(Settings.MaxPrs));
    Channel.W   = Channel.W(   GoodPrs);
    Channel.Prs = Channel.Prs(GoodPrs);

    %find highest and lowest prs level remaining
    z = Channel.Prs(find(Channel.W > 0));
    z = [max(z),min(z)]+[1,-1].*Settings.ZPadding;
  end

  %hence:
  Fine.Prs = single(z(2):Settings.FineGrid(3):z(1));



  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% merge together three 1D functions to make a 3D kernel
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %1d along-track weight
  if Sample.WeightX > 0; Wx = normpdf(Fine.X,0,Sample.WeightX); %simple Gaussian approximation
  else; error('Code not configured for fixed along-track weighting function, stopping'); %not a simple Gaussian - produce from supplied information
  end

  %1d across-track weight
  if Sample.WeightY > 0; Wy = normpdf(Fine.Y,0,Sample.WeightY);  %simple Gaussian approximation
  else;  error('Code not configured for fixed cross-track weighting function, stopping'); %not a simple Gaussian - produce from supplied information
  end

  %1d vertical weight
  if Sample.WeightZ >= 0; z = p2h(10.^Fine.Prs); Wz = normpdf(z,p2h(10.^Sample.Prs),Sample.WeightZ);  %simple Gaussian approximation
  else                                           Wz = interp1(p2h(10.^Channel.Prs),Channel.W,p2h(10.^Fine.Prs),'linear','extrap');%not a simple Gaussian - interpolate from supplied information
  end

  %multiply them together to create a 3D blob
  %order of indexing will seem odd if you're unfamiliar with Matlab
  %don't blame me...
  Wx = repmat(Wx',1,numel(Fine.Y),numel(Fine.Prs));
  Wy = repmat(Wy,numel(Fine.X),1,numel(Fine.Prs));
  Wz = repmat(permute(Wz',[2,3,1]),numel(Fine.X),numel(Fine.Y),1);
  W = Wx.*Wy.*Wz;
  [Fine.Y,Fine.X,Fine.Prs] = meshgrid(Fine.Y,Fine.X,Fine.Prs);

  
  end
