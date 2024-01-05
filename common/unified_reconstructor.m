%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% reconstructor - puts the data into a human-readable shape
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [Error,Output] = unified_reconstructor(Settings,ObsGrid,Final,Simple)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%converts sampled results from list of points to an array in 
%the same format as the original data
%
%Corwin Wright, c.wright@bath.ac.uk, 03/Feb/2023
%
%significantly upated 2023/11/11 to better pass through metadata
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Error = 1; %assume failure unless proved otherwise

%create an output struct
Output = struct();

%for backwards compatability we need to do this:
if isfield(ObsGrid.Recon,'PatternDay'); 
  Output.PatternDay = ObsGrid.Recon.PatternDay;
  ObsGrid.Recon = rmfield(ObsGrid.Recon,'PatternDay'); 
end

%what size should the output be?
DimList = fieldnames(ObsGrid.Recon); 
NDims = numel(DimList); OutputSize = NaN(NDims,1);
for iDim=1:1:NDims; OutputSize(iDim) = max(ObsGrid.Recon.(DimList{iDim}));end; clear iDim

% % %order list of dimensions from largest number of elements to smallest (arbitrary, but I want something stable)
% % %if two have the same size, then order them reverse-alphabetically. This caused me major problems when it occured once unexpectedly...
% % [~,idx] = sortrows(table(OutputSize,fieldnames(ObsGrid.Recon)),[1,2],{'descend' 'descend'});
% % OutputSize = OutputSize(idx);
% % DimList = DimList(idx);
% % clear A;
%actually, just make it alphabetical
[~,idx] = sort(fieldnames(ObsGrid.Recon),'ascend');
OutputSize = OutputSize(idx);
DimList = DimList(idx);
clear A;


%work out where each data point goes
List = ObsGrid.Recon.(DimList{1});
for iDim=2:1:NDims;  List = [List,ObsGrid.Recon.(DimList{iDim})]; end; clear iDim
[~,Order] = sortrows(List,NDims:-1:1);

%drop sampling variables we don't need in the output files
ObsGrid.Track = rmfield(ObsGrid.Track,{'ViewAngleH','ViewAngleZ'});

%reshape all variables with the right number of points into the appropriate ND matrix
%then store them whether we reshaped them or not
Fields = fieldnames(ObsGrid.Track);
for iF=1:1:numel(Fields)
  Var = ObsGrid.Track.(Fields{iF});
  %reshape arrays that have the same number of points as the Recon arrays into the desired shape
  if numel(Var) == numel(Order);ObsGrid.Track.(Fields{iF}) = reshape(Var(Order),OutputSize'); end
  %store output
  Output.(Fields{iF}) = ObsGrid.Track.(Fields{iF});

end
clear Var Fields iF

%also reshape and store our results
Output.T       = reshape( Final(Order),OutputSize');
Output.TSimple = reshape(Simple(Order),OutputSize');

%do we want the output to be singles?
if Settings.SaveSingles == true
  f = fieldnames(Output);
  for iF=1:1:numel(f)
    if strcmpi(class(Output.(f{iF})),'double') %don't convert non-numeric fields
      if ~strcmp(f{iF},'Time') %Matlab times must be stored as doubles
        Output.(f{iF}) = single(Output.(f{iF}));
      end
    end
  end
  clear iF f
end


%done
Error = 0;
return

end
