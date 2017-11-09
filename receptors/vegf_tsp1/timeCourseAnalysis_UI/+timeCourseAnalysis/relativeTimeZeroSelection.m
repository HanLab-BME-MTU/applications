function [ relTimeZero, addTZ, zeroSelect ] = relativeTimeZeroSelection( prompt, default , param)
%% Relative time zero selection
%This is for timecourse analysis
%Obtain relative time zero from user
%User can enter 6 element array [yr month day hr min sec]
% or scalar of the index number of the MD to be used as the relative time zero
% or 'select' to bring up another dialogue box with list dialogue box later
% or 'min' to use MD with earliest acquisition / observation time as time zero

if(nargin < 2 || isempty(default))
%     default = {'2015 7 2 14 44 18.9'};
%     default = {datestr(now,'YYYY mm DD HH MM SS.FFF')};
    default = {'select'};
end

relTimeZero = [];
addTZ = true;
zeroSelect = 0;

userInputStr = inputdlg('6 element time array -or- index number -or- ''select'' -or- ''none''', prompt, 1, default);
if(isempty(userInputStr))
    return;
end
userInputNum = str2num(userInputStr{1});
addTZ = true;

%if time array
if numel(userInputNum) == 6
    relTimeZero = userInputNum;
    zeroSelect = 0; %done no need to do again
end

%if index number
if isscalar(userInputNum)
    if userInputNum > numel(param.fileName)
        error('Input index for relative time zero out of bounds');
    end
    zeroSelect = 1; %need to get relTime0
end

%if not a number
if isempty(userInputNum)
    %if select
    if strcmpi(userInputStr, 'select')
        zeroSelect = 2; %need to get relTime0
    end
    if strncmpi(userInputStr, 'no',2)
        addTZ = false;
        zeroSelect = 0; %done no need to do again
    end
    %if min
    %if strcmpi(userInputStr, 'min')
    %    zeroSelect = 3; %need to get relTime0
    %end
end


% %'same' means same as time above
% userInputStr2 = inputdlg('6 element time array -or- index number -or- ''select'' -or- ''none''', 'Enter VEGF addition time', 1, {'none'});
% if(isempty(userInputStr2))
%     return;
% end
% userInputNum2 = str2num(userInputStr2{1}); %#ok<ST2NM>
% addTZ2 = true;
% %if time array
% if numel(userInputNum2) == 6
%     relTimeZero2 = userInputNum2;
%     zeroSelect2 = 0; %done no need to do again
% end
% %if index number
% if isscalar(userInputNum2)
%     if userInputNum2 > numel(fileName)
%         error('Input index for relative time zero out of bounds');
%     end
%     zeroSelect2 = 1; %need to get relTime0
% end
% %if not a number
% if isempty(userInputNum2)
%     %if select
%     if strcmpi(userInputStr2, 'select')
%         zeroSelect2 = 2; %need to get relTime0
%     end
%     if strcmpi(userInputStr2, 'no VEGF')
%         addTZ2 = false;
%         zeroSelect2 = 0; %done no need to do again
%     end
%     %if min
%     %if strcmpi(userInputStr, 'min')
%     %    zeroSelect = 3; %need to get relTime0
%     %end
% end

%% Relative time zero set
% if index number
if zeroSelect == 1
%     evalc('MD = MovieData.load(movies{userInputNum})');
%     relTimeZero = ML.movies_{userInputNum}.acquisitionDate_;
    zeroSelect = -userInputNum;
end
%if select
if zeroSelect == 2
    zeroSelect = -listdlg('PromptString','Select Movie:', 'SelectionMode','single', 'ListString', param.fileName);
    %evalc('MD = MovieData.load(movies{userChoiceMD})');
%     relTimeZero = ML.movies_{userChoiceMD}.acquisitionDate_;
end

end

