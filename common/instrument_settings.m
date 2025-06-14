function [Settings,InstList] = instrument_settings(Instrument,Settings)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%specific handling options for each observational dataset
%
%contains FineGrid settings for each instrument type, 
%and formatting for calls to load subsets
%
%
%Corwin Wright, c.wright@bath.ac.uk, 01/March/2019
%heavily rewritten 03/February/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%produce a struct containing details of the models, the functions they are called by, and the fine grid settings for sampling that instrument
%'JUSTLISTING' will be an override argument for just getting a list of possible datasets, so don't use this as a name!
%
%%finegrid settings chosen by sensitivity testing as described in WH2018, are ordered as [x,y,Prs], and have units:
%  x: along-track km
%  y: across-track km
%Prs: decade of pressure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%needed for some of the below logic to work when called just to get an instrument list
if nargin== 1; Settings = struct(); Settings.SubSet = ''; end

%create the base struct, then we can start filling it
InstInfo = struct;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set fine grid density for each instrument type
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if ~isnan(sum(Settings.FineGrid)); 

  %override option: set the finegrid to a specific value
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  InstInfo.(Instrument).FineGrid = Settings.FineGrid;

else
  %default finegrids for specified instrument types
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %AIRS (2D)
  InstInfo.AIRS.FineGrid = [2,3,1/20];
  %CRIS (2D)
  InstInfo.CRISj1.FineGrid = [2,3,1/20];  
  InstInfo.CRISsnpp.FineGrid = [2,3,1/20];  
  %AIRS (3D)
  InstInfo.AIRS3D.FineGrid = [1,1,1/20];
  %AIRS (3D) with fake dates
  InstInfo.AIRS3D_fakedates = InstInfo.AIRS3D;
  %COSMIC
  InstInfo.COSMIC.FineGrid = [10,0.5,1/80];
  %HIRDLS
  InstInfo.HIRDLS.FineGrid = [10,2,1/80];
  %limb data, split into regions
  InstInfo.limb_regions.FineGrid = [10,0.5,1/80];
  %gnss resolution test
  InstInfo.gnss_htest = InstInfo.limb_regions;
  %limb, not split into regions
  InstInfo.limb_merged = InstInfo.limb_regions;
  % ALICE/STRIVE
  InstInfo.alice.FineGrid = [20,2,1/80];


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% if we have multiple files per day for an instrument,
% how are the filenames formatted?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Settings.SubSet ~= 0;

  %AIRS (2D)
  InstInfo.AIRS.SubSetInString = ['_g',sprintf('%03d',Settings.SubSet)];

  %CRIS (2D)
  InstInfo.CRISj1.SubSetInString   = ['_g',sprintf('%03d',Settings.SubSet)];
  InstInfo.CRISsnpp.SubSetInString = ['_g',sprintf('%03d',Settings.SubSet)];  
  
  %AIRS (3D)
  InstInfo.AIRS3D.SubSetInString = ['_g',sprintf('%03d',Settings.SubSet)];

  %limb data, split into regions
  InstInfo.limb_regions.SubSetInString = ['_r',sprintf('%03d',Settings.SubSet)];

  % ALICE/STRIVE
  InstInfo.alice.SubSetInString = ['_r',sprintf('%03d',Settings.SubSet)];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%process and return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get list of possible datasets
InstList = fieldnames(InstInfo);

if strcmp(Instrument,'JUSTLISTING')
  %we just want a list of possible instruments - return
  return
else
  % call up the data from the instrument we asked for
  Fields = fieldnames(InstInfo.(Instrument));

  for iF=1:1:numel(Fields)
    Settings.(Fields{iF}) = InstInfo.(Instrument).(Fields{iF});
  end; clear iF Fields

  if ~isfield(Settings,'SubSetInString'); Settings.SubSetInString = ''; end  
  return
end


end

