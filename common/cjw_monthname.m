function MonthID = cjw_monthname(Month,Format)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%takes the number of a month and returns its name
%
%Corwin Wright, corwin.wright@trinity.oxon.org
%13/MAY/2014
%
%inputs
%---------
%
%Month - number (1-12)
%Format - 'full' for full name (default)
%       - 'three' (string) for three-letter identifier
%       - 'one' (string) for one-letter identifier
%
%outputs
%---------
%
%NDays - number of days in month
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2; Format = 'full'; end; 

MonthNames.Full  = {'January','February','March','April','May','June','July','August','September','October','November','December'};
MonthNames.Three = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
MonthNames.One   = {'J','F','M','A','M','J','J','A','S','O','N','D'};

switch Format
  case 'full'
    MonthID = MonthNames.Full{Month};
  case 'three'
    MonthID = MonthNames.Three{Month};
  case 'one'
    MonthID = MonthNames.One{Month};
end

return,MonthID
end
