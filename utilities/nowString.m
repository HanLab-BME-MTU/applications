function str=nowString
%nowString returns a string containing the current date/time without colons or spaces
%
%SYNOPSIS str=nowString
%
%INPUT    none
%
%OUTPUT   str: String containing the current date/time without colons or
%              spaces
%
%c: Jonas, 2/03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

timeStr=datestr(now);
timeStr(12)='-'; %there should be no space in filename
timeStr(findstr(timeStr,':'))='-';

str=timeStr;