function LG_setGoodTag_callback(setGood)
%LG_setGoodTag_callback removes the 2&3 flags from the selected tag or puts
%them back on

% find handles
[naviHandles,movieWindowHandles] = LG_getNaviHandles;

% get tag number
tagIdx = get(gco,'UserData');

% find idlist
idlist = movieWindowHandles.idlist;
goodTimes = movieWindowHandles.idlistData.goodTimes;
dataProperties = movieWindowHandles.dataProperties;

% set good (or bad tag)
[idlist,newDataProperties,dataPropertiesHasChanged] =...
    LG_setGoodTag(idlist,goodTimes,tagIdx,dataProperties,setGood);

% in case dataProperties have changed: Suggest a full recalc
if dataPropertiesHasChanged
    ans = myQuestdlg(['The number of expected spots (maxSpots) has changed.',...
        ' It is strongly suggested to recalculate the idlist'],'Warning',...
        'recalc','reset maxSpots','continue','cancel','recalc');
    switch ans
        case {0,'cancel'}
            % quit
            return
        case 'recalc'
            recalc = 1;
        case 'reset maxSpots'
            recalc = 0;
            newDataProperties = dataProperties;
        case 'continue'
            % don't do anything. Just save
            recalc = 0;
    end
    % set dataProperties
    movieWindowHandles.dataProperties = newDataProperties;
    movieWindowHandles.dataPropertiesHasChanged = dataPropertiesHasChanged;
    guidata(movieWindowHandles.LG_movieWindow,movieWindowHandles);
else
    recalc = 0;
end

if recalc
    idlist = LG_recalcIdlist(idlist,newDataProperties,{0});
end

% save new idlist
LG_loadIdlist(idlist,0);


