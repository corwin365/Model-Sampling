function [Index,Fine,W] = specified_2d_v2(ObsGrid,Instrument,Settings)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%compute fine grid and weight for blobs defined by a specified 2D
%field and a 1D cross-track Gaussian
%
%pre-generate for each unique combination of parameters, and store an
%index, to reduce time in innercore() generating these for each case
%
%
%Corwin Wright, c.wright@bath.ac.uk, 2023/11/11
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%prep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%pull out correct parameters
XZData = ObsGrid.WeightMatrix.XZ.(Instrument);
YSigma = ObsGrid.WeightMatrix.Y.( Instrument);


%choose the right kernel shape for this altitude of tangent point, by interpolating from the surrounding levels
% [~,zidx] = min(abs(p2h(10.^Sample.Prs) - XZData.retrieval_altitude));
% Kernel = XZData.avk(:,:,zidx);
Kernel = interp_1d_ndims(XZData.retrieval_altitude,XZData.avk,p2h(10.^Sample.Prs),3);

x = XZData.distance;
z = XZData.altitude;

%trim weakly-contributing regions
Kernel(abs(Kernel) < (max(abs(Kernel),[],'all') .* Settings.SpecWeightMin)) = 0;
idx = find(sum(Kernel,1) ~= 0); x = x(idx); Kernel = Kernel(:,idx);
idx = find(sum(Kernel,2) ~= 0); z = z(idx); Kernel = Kernel(idx,:);

%find tangent point and set x to be centred here
if isfield(XZData,'tp_altitude')
  %use values provided by Joern (new version)
  [~,idx] = min(abs(XZData.tp_altitude - p2h(10.^Sample.Prs)));
  x = x - XZData.tp_distance(idx);
  clear idx
else
  %estimate TP from total density - old version
  [~,idx] = max(nansum(abs(Kernel),1)); xcentre = x(idx);
  x = x - xcentre;
  clear idx xcentre
end

%fine grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%produce a FineGrid able to represent this field
%X and Y are easy:
Fine.X = min(x):Settings.FineGrid(1):max(x);
Fine.Y = -Settings.BlobScale*YSigma  : Settings.FineGrid(2) : Settings.BlobScale*YSigma ;
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
