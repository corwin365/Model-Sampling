function [Error,Output] = unified_reconstructor(ObsGrid,Final,Simple)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%converts sampled results from list of points to a structure in 
%the same format as the original data
%
%Corwin Wright, c.wright@bath.ac.uk, 03/Feb/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Error = 1; %assume failure unless proved otherwise

%what size should the output be?
DimList = fieldnames(ObsGrid.Recon); NDims = numel(DimList); OutputSize = NaN(NDims,1);
for iDim=1:1:NDims; OutputSize(iDim) = max(ObsGrid.Recon.(DimList{iDim}));end; clear iDim

%order list of dimensions from largest number of elements to smallest (arbitrary, but I want something stable)
[OutputSize,idx] = sort(OutputSize,'desc');DimList = DimList(idx); clear idx

% % % %create output fields
% % % Fields = {'T','TSimple','Lat','Lon','Time','Prs'};
% % % Out = struct();
% % % for iField=1:1:numel(Fields);  Out.(Fields{iField}) = NaN(OutputSize'); end; clear iField

%work out where each data point goes
List = ObsGrid.Recon.(DimList{1});
for iDim=2:1:NDims;  List = [List,ObsGrid.Recon.(DimList{iDim})]; end; clear iDim
[~,Order]  = sortrows(List,NDims:-1:1);

%reshape the metadata
Output = struct();
Fields = {'Lat','Lon','Time','Prs'};
for iField=1:1:numel(Fields)
  Var = ObsGrid.Track.(Fields{iField});
  Var = reshape(Var(Order),OutputSize');
  Output.(Fields{iField}) = Var;
end; clear iField Var

%reshape the output data
Output.T       = reshape( Final(Order),OutputSize');
Output.Tsimple = reshape(Simple(Order),OutputSize');

%done
Error = 0;
return

end
