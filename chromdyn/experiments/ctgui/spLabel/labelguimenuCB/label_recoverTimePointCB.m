function label_recoverTimePointCB
%labelgui menuFunction that allows to recover a timepoint discarded by spotID

%check if there is a movie loaded
imgFigureH = GetUserData(openfig('labelgui','reuse'),'currentWindow');
if isempty(imgFigureH)
    error('no movie loaded.');
    return;
end;

%get time handle
timeH = findall(0,'Tag','slider3');


%try to load data and check whether recovery is possible/allowed
idlist = GetUserData(imgFigureH,'idlist');
slist = GetUserData(imgFigureH,'slist');
if isempty(slist)
    error('sorry, there''s no slist to recover a timepoint from')
    return
end
curr_time = get(timeH,'Value');
dataProperties = GetUserData(imgFigureH,'dataProperties');
idname = GetUserData(imgFigureH,'idname');

if ~isempty(idlist(curr_time).linklist)
    errordlg('idlist is not empty. There''s nothing to recover');
    return
end

if isempty(slist(curr_time).sp)
    errordlg('sorry, timepoint has been rejected by detector already')
    return
end

%we don't want this to happen with idlisttrack and idlisttrack_L
if any(strfind(idname,'idlisttrack'))
    errordlg('sorry, you can''t do that after tracking');
    return
end

if isempty(dataProperties)
    errordlg('can''t recover without dataProperties')
    return
end

%plot data
coord = cat(1,slist(curr_time).sp.cord);
figure(imgFigureH);
plotH = plot(coord(:,1),coord(:,2),'g+');

%find projection and plot there, too
zFigH = GetUserData(imgFigureH,'XZYZFigureH');
if ishandle(zFigH)
    figure(zFigH);
    axesH = GetUserData(zFigH, 'axesH');
    pixelRatio = GetUserData(zFigH, 'pixelRatio');
    movieSize = GetUserData(zFigH, 'movieSize');
    %read and change coordinates (store as coordP, since coord will be used later!)
    coordP = coord;
    coordP(:,3) = (movieSize(3) - coord(:,3)+0.75)*pixelRatio;
    coordP(:,2) = (movieSize(2) - coord(:,2)+1);
    axes(axesH(1));
    plot(coordP(:,1),coordP(:,3),'g+');
    axes(axesH(2));
    plot(coordP(:,2),coordP(:,3),'g+');
end

%ask user whether to recover
ans = questdlg('Do you want to recover this frame (automatic recalc)?','save the spots!','Yes','No','Yes');

if isempty(ans)
    delete(plotH);
    return
    %user cancelled
end

switch ans
    case 'No'
        %remove plotted 
        delete(plotH);
        return
    case 'Yes'
        %fill recoverSlist with info
        PIXELSIZE_XY=dataProperties.PIXELSIZE_XY;
        PIXELSIZE_Z=dataProperties.PIXELSIZE_Z;
        recoverSlist.nspots=size(coord,1);
        if recoverSlist.nspots>10
            error('more than 10 spots. Can''t be handles by spotID')
            return
        end
        p2m=ones(recoverSlist.nspots,1)*[PIXELSIZE_XY PIXELSIZE_XY PIXELSIZE_Z];
        recoverSlist.xyz = coord.*p2m;
        recoverSlist.amp = cat(1,slist(curr_time).sp.amp);
        if ~isempty(slist(curr_time).statistics)
            recoverSlist.detectQ=slist(curr_time).statistics.Q;
            recoverSlist.noise=slist(curr_time).statistics.chi;
        else
            recoverSlist.detectQ=[];
            recoverSlist.noise=[];
        end
        recoverSlist.trackQ=[]; 
        recoverSlist.time = curr_time;
        recoverSlist.CoM = mean(recoverSlist.xyz,1);
end


%idlist needs to be recalced from the beginning
idlist = recalcIdlist(idlist,1,[],dataProperties,recoverSlist);

%save data
SetUserData(imgFigureH,idlist,1);

%get view3DH if exist and update data
view3DH = GetUserData(imgFigureH,'view3DGUIH');
if ishandle(view3DH)
    view3D_generateHandles;
end

%refresh labelgui
labelgui('refresh');
