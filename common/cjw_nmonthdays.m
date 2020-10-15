function [NDays] = cjw_nmonthdays(Month,Year)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%give number of days in a given month
%takes several month formats; year optional (assumes non-leap year)
%
%Corwin Wright, corwin.wright@trinity.oxon.org
%13/JAN/2014
%
%inputs
%---------
%
%Month - can be number, three-letter ID, or full name (case irrelevant)
%Year  - four-digit format, optional
%
%outputs
%---------
%
%NDays - numnber of days in month
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2; Year = 2013; end; %arbitrary non-leap year


%if we fed in a string, try and identify the month from it
if strcmp(class(Month),'char') == 1;
  
  Month = lower(Month); %convert to lower case
  
  switch Month
    case 'jan';      Month = 1;
    case 'feb';      Month = 2;
    case 'mar';      Month = 3;
    case 'apr';      Month = 4; 
    case 'may';      Month = 5;
    case 'jun';      Month = 6;
    case 'jul';      Month = 7;
    case 'aug';      Month = 8; 
    case 'sep';      Month = 9;
    case 'oct';      Month = 10;
    case 'nov';      Month = 11;
    case 'dec';      Month = 12;
    case 'january';  Month = 1;
    case 'february'; Month = 2;
    case 'march';    Month = 3;
    case 'april';    Month = 4; 
    %may is the same in three-letter or full name!
    case 'june';     Month = 6;
    case 'july';     Month = 7;
    case 'august';   Month = 8; 
    case 'september';Month = 9;
    case 'october';  Month = 10;
    case 'november'; Month = 11;
    case 'december'; Month = 12;        
    otherwise;       Month = NaN;
  end
  
end

%leap year?
%may as well cover all our bases here...
if     mod(Year,400) == 0; LeapYear = 1;
elseif mod(Year,100) == 0; LeapYear = 0;
elseif mod(Year,  4) == 0; LeapYear = 1;
else                       LeapYear = 0;
end


%error checking
if isnan(Month); NDays = 0;
else
  MonthLengths = [31,28,31,30,31,30,31,31,30,31,30,31];
  NDays = MonthLengths(Month);
  if LeapYear == 1; if Month == 2; NDays=NDays+1; end; end;
end





