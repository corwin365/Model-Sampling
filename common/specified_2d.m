  function [Fine,W] = specified_2d(Sample,ObsGrid,Instrument,Settings)


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %compute fine grid and weight for blobs defined by a specified 2D
  %field and a 1D cross-track Gaussian
  %
  %Corwin Wright, c.wright@bath.ac.uk, 2023/11/11
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %prep
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %pull out correct parameters
  XZData = ObsGrid.WeightMatrix.XZ.(Instrument);
  YSigma  = ObsGrid.WeightMatrix.Y.( Instrument);

  % % % %interpolate betwen levels to find the right kernel for this level
  % % % %if you don't do this then the results show steps between heights
  % % % zo = p2h(10.^Sample.Prs);
  % % % zi = XZData.retrieval_altitude;
  % % % Kernel = interp_1d_ndims(zi,XZData.avk,zo,3);

  %choose the right kernel shape for this altitude of tangent point
  [~,zidx] = min(abs(p2h(10.^Sample.Prs) - XZData.retrieval_altitude));
  Kernel = XZData.avk(:,:,zidx);
  x = XZData.distance;
  z = XZData.altitude;

  %trim weakly-contributing regions
  Kernel(abs(Kernel) < max(abs(Kernel),[],'all') .* Settings.SpecWeightMin) = 0;
  idx = find(sum(Kernel,1) ~= 0); x = x(idx); Kernel = Kernel(:,idx);
  idx = find(sum(Kernel,2) ~= 0); z = z(idx); Kernel = Kernel(idx,:);  

  %temporarily adjust kernel to be centred at TP - Joern's next version will fix this
  [~,idx] = max(nansum(abs(Kernel),1)); xcentre = x(idx);
  x = x - xcentre; 

  %fine grid
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %produce a FineGrid able to represent this field
  %X and Y are easy:
  Fine.X = min(x):Settings.FineGrid(1):max(x);
  Fine.Y = -1.25*Settings.BlobScale*YSigma  : Settings.FineGrid(2) : 1.25*Settings.BlobScale*YSigma ;
  Fine.Prs = log10(h2p(max(z))):Settings.FineGrid(3):log10(h2p(min(z)));

  %interpolate XZ kernel onto the finegrid
  z = log10(h2p(z));
  [xo,yo] = meshgrid(Fine.X,Fine.Prs);
  [xi,yi] = meshgrid(x,z);
  Kernel = interp2(xi,yi,Kernel,xo,yo);


 
  %weight
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %1d across-track weight
  Wy = normpdf(Fine.Y,0,YSigma);  %simple Gaussian approximation

  %multiply with the XZ field to to create a 3D blob
  Kernel = repmat(permute(Kernel,[1,3,2]),[1,numel(Fine.Y),1]);
  Wy     = permute(repmat(Wy,[numel(Fine.X),1,numel(Fine.Prs)]),[3,2,1]);
  W      = permute(Kernel.*Wy,[3,2,1]);

  %done!
  [Fine.Y,Fine.X,Fine.Prs] = meshgrid(Fine.Y,Fine.X,Fine.Prs);


 end