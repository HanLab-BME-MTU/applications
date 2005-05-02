function label_removeFusionsCB
%LABEL_REMOVEFUSIONSCB is a callback to remove all fusion tags in a movie
%


% collect data

cffig = openfig('labelgui','reuse');

imgFigureH = GetUserData(openfig('labelgui','reuse'),'currentWindow');
if isempty(imgFigureH)
    return;
end;
dataProperties = GetUserData(imgFigureH,'dataProperties');
idlist = GetUserData(imgFigureH,'idlist');

% find frames to remove: First collect the number of spots in every frame.
% Wherever we have less than the maximum number of spots, there is a
% fusion.
nTimepoints = length(idlist);
% fill spotlist
nSpots = zeros(nTimepoints,1);
for t=1:nTimepoints
    if ~isempty(idlist(t).linklist)
        nSpots(t) = max(idlist(t).linklist(:,2));
    end
end

% index of spots to remove - there's nothing to remove if there isn't
% anything.
nTags = max(nSpots);
rmIdx = find(nSpots < nTags & nSpots > 0);
        

ButtonName = questdlg(sprintf('You are about to remove %i/%i frames',length(rmIdx),length(find(nSpots))),'WARNING','Remove','Cancel','Cancel');
if strcmp(ButtonName,'Remove')
    %continue
else
    return
end

% remove spots. 
[idlist(rmIdx).linklist] = deal([]);

% relink spots
% the easiest way for this is to just rerun spotID - without all the
% fusions, we shouldn't have a problem
idlist = recalcIdlist(idlist,1,[],dataProperties);

% store idlist and refresh
%save data
SetUserData(imgFigureH,idlist,1);

%get view3DH if exist and update data
view3DH = GetUserData(imgFigureH,'view3DGUIH');
if ishandle(view3DH)
    view3D_generateHandles;
end

%update labelgui
labelgui('refresh');
