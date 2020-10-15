function [FakeData,OldFile] = extract_fake_data(MatlabDay,FileString,OldFile)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%extract 'fake' (i.e. reanalysis or model sampled as satellite) data 
%for a given day using the standard format for this project
%only runs if calling a new file to previous iteration
%
%Corwin Wright, c.wright@bath.ac.uk
%09/NOV/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%errors:
%0. no error - successful
%1. no file found
%2. no data on this day
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find the file (if it's not the one from last time) and load it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(OldFile,'Name'); OldFile.Name = ' '; end; %always load file afresh if old file not specified


if exist(FileString) == 0;   
  FakeData.Error = 1;
  return;
end

if strcmp(FileString,OldFile.Name) == 0;
  %new file - load it up
   AllData = load(FileString);
   AllData = AllData.Sampled_Data;

 
   %store the name, to avoid need to reload
   OldFile.Name = FileString;
  
else
  %same as last call - don't reload
  AllData = OldFile.Data;
end
clear m y DataDir FileString

OldFile.Data = AllData;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%extract the day we actually want
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if numel(AllData.T) ==0; 
  FakeData.Error = 2;
  return;
end
clear MatlabDay

FakeData.Data = AllData;
FakeData.Error = 0; %it worked!

return
end
