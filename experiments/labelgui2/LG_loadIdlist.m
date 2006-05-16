function success = LG_loadIdlist(idlist, replace, idname, idlistDir, dataFileName)
%LG_loadIdlist writes a new idlist into the labelgui
%
% replace = 1 means that we have to discard the old safeIdlist
% idlistDir is optional if replace = 0
% LG_loadIdlst will always plot at the end.


[naviHandles,movieWindowHandles] = LG_getNaviHandles;

% read data from idlist. IdlistData has the fields:
% nSpots;
% goodTimes
% maxSpots
% maxTags
% labelcolor
% labellist
% flagList
idlistData = LG_readIdlistData(idlist,movieWindowHandles.dataProperties);

if replace
    % now that we know the maximum number of tags, we can set a colormap. The
    % first maxNumSpots tags get distributed evenly, then we arrange the rest.
    % Do not update colorMap if we're just updating the idlist
    cmapLength = idlistData.maxTags + ~isEven(idlistData.maxTags);
    firstIndices = [1:floor(...
        cmapLength/...
        (movieWindowHandles.dataProperties.MAXSPOTS+1)):cmapLength];
    otherIndices = missingIndices(firstIndices,cmapLength)';
    cmap = hsv(cmapLength);
    % rearrange
    cmap = cmap([firstIndices,otherIndices(1:2:end),otherIndices(2:2:end)],:);
    movieWindowHandles.colorMap = cmap;
else
    % if we're updating the idlist, it is possible that we recalc'ed, and
    % thereby increased the number of entries in linklist. In this case, we
    % will just reuse the first unused color of cmap for now
    cmap = movieWindowHandles.colorMap;
    nColors = size(cmap,1);
    if idlistData.maxTags > nColors
        if idlistData.maxSpots < nColors
            % reuse the last color in cmap. Remember maxSpots entries
            cmapNew = repmat(cmap(end,:),idlistData.maxTags,1);
            cmapNew(1:idlistData.maxSpots,:) = ...
                cmap(1:idlistData.maxSpots,:);
        elseif idlistData.maxSpots == nColors
            % Use 001 for new tags
            cmapNew = repmat([0,0,1],idlistData.maxTags,1);
            cmapNew(1:nColors,:) = ...
                cmap;
        else
            % I don't think this should happen. If it does, give warning
            % and totally reassing cmap

            warning('unexpected reassignment of colormap in LG_loadIdlist. Please inform support')

            cmapLength = idlistData.maxTags + ~isEven(idlistData.maxTags);
            firstIndices = [1:floor(...
                cmapLength/...
                (movieWindowHandles.dataProperties.MAXSPOTS+1)):cmapLength];
            otherIndices = missingIndices(firstIndices,cmapLength)';
            cmap = hsv(cmapLength);
            % rearrange
            cmapNew = cmap([firstIndices,otherIndices(1:2:end),otherIndices(2:2:end)],:);

        end
        movieWindowHandles.colorMap = cmapNew;
    end
end


% store idlist, dataFileName, idlistDir, set dataHasChanged to 0
% also store a copy of the idlist in case we need to revert
movieWindowHandles.idlist = idlist;
if replace
    movieWindowHandles.dataHasChanged = 0;
    movieWindowHandles.safeIdlist = idlist;
    % idlistDir is the directory we loaded the idlist from. Normally, it is the
    % same as movieDir, of course
    movieWindowHandles.idlistDir = idlistDir;
    movieWindowHandles.dataFileName = dataFileName;
    movieWindowHandles.idname = idname;
else % if we aren't replacing, idlist will be different from safeIdlist
    movieWindowHandles.dataHasChanged = 1;
end
movieWindowHandles.idlistData = idlistData;


% write context menu
movieWindowHandles = LG_createTagPopupMenu(movieWindowHandles,...
    movieWindowHandles.LG_movieWindow);


% Flag. If replace, we set to all frames. Otherwise, keep current flag
% selection and call the callback to update flaggedFrameList
flagNames = LG_getFlagNames(idlistData.flagList);
movieWindowHandles.flagData.flagNames = flagNames;

if replace
    flaggedFrameList = (1:movieWindowHandles.dataProperties.movieSize(4))';
    flagIdx = LG_getCurrentTime;

    % set flag-PD
    set(naviHandles.LG_navi_flagName_pd,'String',flagNames(:,1),'Value',1);

    % store flagNames, flagFrameList
    movieWindowHandles.flagData.flaggedFrameList = flaggedFrameList;
    movieWindowHandles.flagData.flagIdx = flagIdx;

    % and save everything
    guidata(movieWindowHandles.LG_movieWindow,movieWindowHandles);

else

    % keep old value and update via callback - save everything first,
    % though
    guidata(movieWindowHandles.LG_movieWindow,movieWindowHandles);

    % we don't want the function to plot
    LG_navi_flagName_pd_Callback(0);

end


% now that all is loaded, we can go and plot current frame.
if replace
    LG_gotoFrame(1);
else
    LG_gotoFrame(LG_getCurrentTime);
end

% update reAssignGUI
if ~isempty(movieWindowHandles.otherWindows.LG_reAssignGUI)
    close(movieWindowHandles.otherWindows.LG_reAssignGUI);
    LG_reAssignGUI;
end

% launch all the view-windows
if strcmp(get(naviHandles.LG_navi_menuShowTestRatios,'checked'),'on');
    set(naviHandles.LG_navi_menuShowTestRatios,'checked','off');
    LG_showTestRatios_callback;
end
if strcmp(get(naviHandles.LG_navi_menuShowIntensities,'checked'),'on');
    if replace      
    set(naviHandles.LG_navi_menuShowIntensities,'checked','off');
    LG_showIntensities_callback;
    else
        LG_showIntensities_callback(1);
    end
end
if strcmp(get(naviHandles.LG_navi_menuShowDistances,'checked'),'on');
    set(naviHandles.LG_navi_menuShowDistances,'checked','off');
    LG_showDistances_callback;
end
if strcmp(get(naviHandles.LG_navi_menuShowDisplacements,'checked'),'on');
    set(naviHandles.LG_navi_menuShowDisplacements,'checked','off');
    LG_showDisplacements_callback;
end

success = 1;