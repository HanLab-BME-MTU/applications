function LG_navi_flagName_pd_Callback(doPlot)
%LG_navi_flagName_pd_Callback is the callback for the flag-pulldown

if nargin == 0 || isempty(doPlot)
    doPlot = 1;
end

% get handles and data
[naviHandles,movieWindowHandles] = LG_getNaviHandles;
flagList = movieWindowHandles.idlistData.flagList;
flagNames = movieWindowHandles.flagData.flagNames;

pdH = naviHandles.LG_navi_flagName_pd;

% find the chosen flag
currentFlag = flagNames{get(pdH,'Value'),2};

% switch according to flags. Build new flaggedFrameList, set idx
switch currentFlag
    case -2 % non-flagged frames
        flaggedFrameList = find(all(flagList==0,2));
    case -1 % all frames
        flaggedFrameList = [1:size(flagList,1)];
    case 0 % all flagged frames (no deleted frames)
        flaggedFrameList = find(any(flagList > 0,2) & all(flagList < 21,2));
    case {1,2,3,11,12,13,21} % particular flags
        flaggedFrameList = find(any(flagList == currentFlag,2));
    otherwise
        h = errordlg('Flag option not considered','Incomplete code');
        uiwait(h);
        return
end

% store data
movieWindowHandles.flagData.flaggedFrameList = flaggedFrameList;
guidata(movieWindowHandles.LG_movieWindow, movieWindowHandles);


% if we plan to plot, let gotoFlag take care of updating the flagNumber.
% Else, update flag number
if doPlot
    % plot.
    LG_gotoFlag(0);
else
    if isempty(flaggedFrameList)

        % remove flagNumber
        set(naviHandles.LG_navi_flagNumber_txt,...
            'String','');


    else
        % write flagNumber
        set(naviHandles.LG_navi_flagNumber_txt,...
            'String',sprintf('%i/%i',1,length(flaggedFrameList)));

    end
end

