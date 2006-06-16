function LG_navi_flagName_pd_Callback(doPlot)
%LG_navi_flagName_pd_Callback is the callback for the flag-pulldown

% {'All Frames',         -1;...
%     'Flagged Frames',       0;...
%     'Estimated Spots',      1;...
%     'Tracked Spots',        2;...
%     'Fusion Spots',         3;...
%     'Single Occurences',   12;...
%     'Deleted Frames',      21;...
%     'Tags Close To Border',22;...
%     'Adjusted Intensities',11;...   
%     'Non-flagged Frames',  -2;...
%     };

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
        flaggedFrameList = 1:size(flagList,1);
    case 0 % all flagged frames [tracked, fusions, single occ., close to border]
        flaggedFrameList = find(any(flagList == 2 | flagList == 3 | flagList == 12 | flagList == 22, 2))
    case 2 % all tracked frames
        flaggedFrameList = find(any(flagList == 2 | flagList == 5,2));
    case 3 % all fusion frames
        % primary fusions can be 3 or 5; secondary fusions are always there
        % too, and always 4
        flaggedFrameList = find(any(flagList == 4, 2));
    case {1,11,12,13,21,22} % particular flags
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

