function LG_gotoFlag(direction)
%LG_gotoFlag sets currentTime to previous, current or next flag

% I originally wanted to store the index of the flag that was last visited,
% so that after going a few frames ahead and a few back, going to the next
% flag would not bring you to the one you just visited. However, this is
% problematic if you, say, delete the current, flagged frame, because then
% the previous flag#2 will become flag#1 and advancing by a flag will bring
% you to what was flag#3. Therefore, I decided to just go to the next
% possible flagged frame

% input: direction
% 1 go to closest flag later than current time
% 0 go to first flag
% -1 go to closest flag earlier than current time

% flagNames = ...
%     {'All Frames',         -1;...
%     'Flagged Frames',       0;...
%     'Estimated Spots',      1;...
%     'Tracked Spots',        2;...
%     'Fusion Spots',         3;...
%     'Single Occurences',   12;...
%     'Deleted Frames',      21;...
%     'Adjusted Intensities',11;...
%     'Non-flagged Frames',  -2;...
%     };

% get handles
[naviHandles, movieWindowHandles] = LG_getNaviHandles;

if isempty(movieWindowHandles)
    return
end

% get flag data
flaggedFrameList = movieWindowHandles.flagData.flaggedFrameList;

nFlags = length(flaggedFrameList);

% get current time
currentTime = LG_getCurrentTime;


if nFlags == 0
    % remove flagNumber
    set(naviHandles.LG_navi_flagNumber_txt,...
        'String','');
    % turn labelgui invisible so that the user first closes the dialogue
    % before continuing
    set(naviHandles.LG_navigator,'Visible','off');
    % inform user
    h = warndlg(sprintf('There is no frame with this flag!'));
    uiwait(h)
    set(naviHandles.LG_navigator,'Visible','on');
    return
end

% Move and plot
switch direction
    case +1
        % go forward one flag
        delta = flaggedFrameList - currentTime;
        % use ML7 feature with find. The frameList is ordered, so the next
        % flag will be the first where flagTime>currentTime
        nextIdx = find(delta > 0,1);

        % update only if possible!
        if ~isempty(nextIdx)
            % plot
            LG_gotoFrame(flaggedFrameList(nextIdx));
        end


    case -1
        % go back one flag
        delta = flaggedFrameList - currentTime;
        % use ML7 feature with find. The frameList is ordered, so the next
        % flag will be the last where flagTime<currentTime
        prevIdx = find(delta < 0,1,'last');

        % update only if possible!
        if ~isempty(prevIdx)
            % plot
            LG_gotoFrame(flaggedFrameList(prevIdx));
        end


    case 0
        % goto first flag
        LG_gotoFrame(flaggedFrameList(1));
end
