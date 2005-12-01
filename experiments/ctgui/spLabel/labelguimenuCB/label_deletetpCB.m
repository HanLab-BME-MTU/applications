function label_deletetpCB

cffig = openfig('labelgui','reuse');

imgFigureH = GetUserData(openfig('labelgui','reuse'),'currentWindow');
if isempty(imgFigureH)
    return;
end;
dataProperties = GetUserData(imgFigureH,'dataProperties');

ButtonName = questdlg('Remove current time point?','Warning','Yes','Yes&Recalc','No','No');
switch strcmp(ButtonName,'No')+2*(strcmp(ButtonName,'Yes'))
case 1 %no
    return %end execution here
case 2 %yes
    recalc = 0;
case 0 %yes%recalc
    recalc = 1;
end

%read data
timeslideH = findall(cffig,'Tag','slider3');
idlist = GetUserData(imgFigureH,'idlist');
curr_time = get(timeslideH,'Value');

% remove curr_time
idlist(curr_time).linklist = [];

%find prev_time
done = 0;
step = 1;
prev_time = [];
while ~done&(curr_time-step>0)
    if ~isempty(idlist(curr_time-step).linklist)
        prev_time = curr_time-step;
        done = 1;
    end
    step = step+1;
end

%find next time
done = 0;
step = 1;
next_time = [];
while ~done&(curr_time+step<length(idlist))
    if ~isempty(idlist(curr_time+step).linklist)
        next_time = curr_time+step;
        done = 1;
    end
    step = step+1;
end

%remap linklists
switch (isempty(prev_time))+2*(isempty(next_time))
case 0 %both exist
    for i = 1:size(idlist(prev_time).linklist,1)
        idlist(prev_time).linklist(i,7) = ...
            idlist(next_time).linklist(find(idlist(prev_time).linklist(i,4)==idlist(next_time).linklist(:,4)),2);
        idlist(next_time).linklist(i,6) = ...
            idlist(prev_time).linklist(find(idlist(next_time).linklist(i,4)==idlist(prev_time).linklist(:,4)),2);
    end
    tStart = prev_time;
    
case 1 %no prev_time
    idlist(next_time).linklist(:,6) = 0;
    tStart = curr_time;
    
case 2 %no next_time
    idlist(prev_time).linklist(:,7) = 0;
    tStart = prev_time; %actually, no linking will be done
end

%write idlist-status
idlist(1).stats.status{end+1}=[date,': deleted frame ',num2str(curr_time)];


%recalculate connections
if recalc
    idlist = recalcIdlist(idlist,tStart,[],dataProperties);
end
        
%save data
SetUserData(imgFigureH,idlist,1);

%get view3DH if exist and update data
view3DH = GetUserData(imgFigureH,'view3DGUIH');
if ishandle(view3DH)
    view3D_generateHandles;
end

%update labelgui
labelgui('refresh');

