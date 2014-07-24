function label_recalcconnections
%menu item for labelgui: recalculate connections between timepoint with reCalCon
%12/02/JD

%get handles
imgFigureH = GetUserData(openfig('labelgui','reuse'),'currentWindow');
if isempty(imgFigureH)
    return;
end;
timeH = findall(0,'Tag','slider3');


%load data
idlist = GetUserData(imgFigureH,'idlist');
curr_time = get(timeH,'Value');
dataProperties = GetUserData(imgFigureH,'dataProperties');


prompt   =  {'Start at frame','End at frame','Weight between Start and End (Values less than .5 favor intensity Values over .5 favor displacement)'};
title    =  'Recalculate connections';
lines =  1;
def      =  {num2str(curr_time),num2str(size(idlist,2)),num2str(idlist(1).stats.weight(1))};
userInput   =  inputdlg(prompt,title,lines,def);


%cancel?
if isempty(userInput)
    return %end evaluation here
end

tstart = str2num(char(userInput{1}));
tend = str2num(char(userInput{2}));
opt.weight = str2num(char(userInput{3}));

%correct weight?
 if opt.weight<0 | opt.weight>1
     h = warndlg('Wrong weight (has to be 0...1)','Warning');
     uiwait(h);
     label_recalcconnections;
     return
 end

%recalc
idlist = recalcIdlist(idlist,[tstart,tend],opt,dataProperties);

%save data
SetUserData(imgFigureH,idlist,1);

%get view3DH if exist and update data
view3DH = GetUserData(imgFigureH,'view3DGUIH');
if ishandle(view3DH)
    view3D_generateHandles;
end

%refresh labelgui
labelgui('refresh');



