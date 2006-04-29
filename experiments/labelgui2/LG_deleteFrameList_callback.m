function LG_deleteFrameList_callback
%LG_deleteFrame is the callback for deleting a list of frames

% get handles
[naviHandles, movieWindowHandles] = LG_getNaviHandles;

% get idlist
idlist = movieWindowHandles.idlist;
dataProperties = movieWindowHandles.dataProperties;

% find currentTime, goodTimes
currentTime = LG_getCurrentTime;
goodTimes = movieWindowHandles.idlistData.goodTimes;

% ask user for list of frames to remove (first, last)
% also allow to set recalc (y,n). Add maybe more options in the future,
% e.g. rm all frames with fusions etc.
prompt = {'Start frame','End frame','Recalc (0=No/1=Yes)?'};
name = 'remove frame list';
defaultAnswer = {num2str(currentTime),num2str(length(idlist)),'0'};
answer = inputdlg(prompt,name,1,defaultAnswer);

if isempty(answer) 
    return
end

% read answer
firstFrame = str2double(answer{1});
lastFrame = str2double(answer{2});
recalc = str2double(answer{3});

% make listOfFrames
listOfFrames = firstFrame:lastFrame;
if isempty(listOfFrames)
    return
end


% delete frame
idlist = LG_deleteFrame(idlist,dataProperties,listOfFrames,goodTimes,recalc);


% save idlist, idlistData - simply reload without replacing everything
LG_loadIdlist(idlist, 0);
