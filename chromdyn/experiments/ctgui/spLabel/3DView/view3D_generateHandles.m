function view3D_generateHandles
%generating handle structure for links and points


%get necessary handles and data
labelPanelH = GetUserData(openfig('labelgui','reuse'),'currentWindow');
idlist = GetUserData(labelPanelH,'idlist');
labelguiH = findall(0,'Tag','labelgui');
view3DH = findall(0,'Tag','view3DGUI');
view3D_handles = guidata(view3DH); %handle structure of view3D

cMap = GetUserData(labelguiH,'cMap');
cMapFact = size(cMap,1)/idlist(1).stats.maxColor;
tmax = size(idlist,2);

%init myWaitbar, use 2*tmax because 2 loops are used
maxval = 2*tmax;
waitbarHandle = mywaitbar(0,[],maxval,'building list of object handles...');

%------------------build spot list, taking into account centering/alignment----------



%find first good timepoint
t1 = 1;
if isempty(idlist(t1).linklist)
    %lookfor first good frame
    done = 0;
    while ~done
        t1 = t1+1;
        if ~isempty(idlist(t1).linklist);
            done = 1;
        end
    end
end


%coordinates switchboard

%centering
centerMode = get(view3D_handles.view3D_center_PD,'Value');
switch centerMode
    case 1 %center none
        centerString = '';
    case 2 %center centroid
        centerString = '-idlist(t).centroid';
    case 3 %center tag
        tagNum = view3D_handles.centerTagNumber; % = log2(tagColor)+1
        centerString = ['-idlist(t).linklist(',num2str(tagNum),',9:11)'];
end

%aligning
alignMode = get(view3D_handles.view3D_align_PD,'Value');
switch alignMode
    case 1 %align none
        alignString = '';
        transposeString = '';
    case 2 %align axis
        rotMatrix = view3D_handles.rotMatrix;
        alignString = 'rotMatrix(:,:,t)*';
        transposeString = '''';
    case 3 %align metaphase
        rotMatrixMeta = view3D_handles.rotMatrixMeta;
        alignString = 'rotMatrixMeta(:,:,t)*';
        transposeString = '''';
end



%build coord-assignment string
coordString = ['slist(t).spotlist(i,4:6) = ',alignString,'(idlist(t).linklist(rowIdx(1),9:11)',centerString,')',transposeString,'; '];

nspots = zeros(tmax,1);
slist(1:tmax) = struct('spotlist',[],'linkdown',[]);


try
    for t = 1:tmax
        if ~isempty(idlist(t).linklist)
            %sort idlist (does not need separate loop)
            idlist(t).linklist = sortrows(idlist(t).linklist,4);
            nspots(t) = max(idlist(t).linklist(:,2));
            for i = 1:nspots(t);
                rowIdx = find(idlist(t).linklist(:,2)==i);
                slist(t).spotlist(i,1) = t;
                slist(t).spotlist(i,2) = i;
                eval(coordString); %coordinates
                slist(t).spotlist(i,3) = (idlist(t).linklist(rowIdx(1),3)/length(rowIdx)); %color, already adjusted for correct fusion
                slist(t).linkdown{i} = unique(idlist(t).linklist(rowIdx,7))'; %linkdown
                slist(t).spotlist(i,7) = sum(idlist(t).linklist(rowIdx,8)); %amplitude
            end
        end
        mywaitbar(t/maxval,waitbarHandle,maxval);
    end
catch
    if findstr(lasterr,['Error using ==> get',char(10),'Invalid handle'])
        error('evaluation aborted by user');
    else
        rethrow(lasterror) 
    end
end

goodTime = find(nspots);

%calculate plot size (extreme coordinates)
allSpots = cat(1,slist.spotlist);
xCoord(2,:) = max(allSpots(:,4:6),[],1);
xCoord(1,:) = min(allSpots(:,4:6),[],1);
centerC = mean(xCoord,1); %center coordinate
xCoordAir = (xCoord-[centerC;centerC])*1.05+[centerC;centerC];%add 5% 'air' around plot: [0;1]->[-.5;.5]->[-.5025,.5025]->[-.025,1.025]

%if aligning, it can happen that all motion is confined to a plane. This
%does not look good; neither do the plots with very unequal axis lenghts.
%therefore: Axis-limits are set to min(xCoordAir(1,:)/max(xCoordAir(2,:)))

minAxLimit = min(xCoordAir(1,:));
maxAxLimit = max(xCoordAir(2,:));


%calculate min, max intensity
xInt(2) = max(allSpots(:,7),[],1);
xInt(1) = min(allSpots(:,7),[],1);

%calculate parameters to adjust tag marker size: smallest intensity->1/30
%of axis, highest intensity->1/15 of axis
shAxis = (maxAxLimit-minAxLimit)/30;
deltaXInt = xInt(2)-xInt(1);

%create figure if there isn't one already
view3DfigH = GetUserData(labelPanelH,'view3DfigH');
if isempty(view3DfigH)
    view3DfigH = figure;
    set(view3DfigH,'Tag','view3Dfig','Name','3DViewFigure','HandleVisibility','Callback');
end
%save figure handle
SetUserData(labelPanelH,view3DfigH,1);


%initialize figure
figure(view3DfigH);
cla; %clear previous data
axesH = gca;
set(axesH,'NextPlot','add');
set(axesH,'XGrid','on');
set(axesH,'YGrid','on');
set(axesH,'ZGrid','on');
set(axesH,'XLim',[minAxLimit,maxAxLimit]);
set(axesH,'YLim',[minAxLimit,maxAxLimit]);
set(axesH,'ZLim',[minAxLimit,maxAxLimit]);
daspect([1,1,1]);;
plot3(centerC(1),centerC(2),centerC(3),'w');
xlabel('x [\mum]');
ylabel('y [\mum]');
zlabel('z [\mum]');

%plot invisible graphics&store handles in cell arrays
% spotData{1,1} = [];
% spotDataTmp = [];
arrowHandles = cell(tmax,1);
colorList = cell(tmax,1);
arrowHandlesTmp = [];
colorListTmp = [];

%turn of warning divide by zero to be able to plot 0-arrows
oldWarnings = warning;
warning off MATLAB:divideByZero;

try
    for t2 = t1+1:tmax
        figure(view3DfigH);
        
        if ~isempty(slist(t2).spotlist)
            for i = 1:nspots(t1)
                %spot color
                spotCol = cMap(round(slist(t1).spotlist(i,3)*cMapFact),:);
                %draw all arrows starting from spot(t1,i)
                linkdownList = slist(t1).linkdown{i};
                for j = 1:length(linkdownList)
                    %get starting and end points, then calculate vectorlength to determine size of arrowhead
                    %used right now: head-height:width = 2:1, height = 0.2*length
                    arrowStartEnd = [slist(t1).spotlist(i,4:6);slist(t2).spotlist(linkdownList(j),4:6)];
                    arrowLength = norm(arrowStartEnd(2,:)-arrowStartEnd(1,:));
                    arrowH = arrow3(arrowStartEnd(1,:),arrowStartEnd(2,:),...
                        [],.1*arrowLength,.2*arrowLength,((slist(t1).spotlist(i,7)-xInt(1))/deltaXInt+1)*shAxis);
                    set(arrowH(1),'Color',spotCol);
                    set(arrowH(2:3),'FaceColor',spotCol);
                    %set(arrowH,'Visible','off');
                    arrowHandlesTmp = [arrowHandlesTmp;arrowH'];
                    colorListTmp = [colorListTmp;spotCol];
                end
            end
            arrowHandles{t1} = arrowHandlesTmp;
            arrowHandlesTmp = [];
            colorList{t1} = colorListTmp;
            colorListTmp = [];
            t1 = t2;
        end
        mywaitbar((t2+tmax)/maxval,waitbarHandle,maxval);
    end
    %for last timepoint: plot only spot
    for i = 1:nspots(t1)
        %spot color
        spotCol = cMap(round(slist(t1).spotlist(i,3)*cMapFact),:);
        %draw all arrows starting from spot(t1,i)
        linkdownList = slist(t1).linkdown{i};
        for j = 1:length(linkdownList)
            arrowH = arrow3(slist(t1).spotlist(i,4:6),slist(t1).spotlist(i,4:6),...
                [],.1,.1,slist(t1).spotlist(i,7)*0.01);
            set(arrowH(3),'FaceColor',spotCol);
            %set(arrowH,'Visible','off');
            arrowHandlesTmp = [arrowHandlesTmp;arrowH'];
            colorListTmp = [colorListTmp;spotCol];
        end
    end
catch
    if findstr(lasterr,['Error using ==> get',char(10),'Invalid handle'])
        error('evaluation aborted by user');
    else
        rethrow(lasterror) 
    end
end

arrowHandles{t1} = arrowHandlesTmp;
arrowHandlesTmp = [];
colorList{t1} = colorListTmp;
colorListTmp = [];
%turn on old warnings
warning(oldWarnings);

%if there is axisdata: update (save current handles first!)
if view3D_handles.axisDefined~=0
    guidata(view3D_handles.view3DGUI,view3D_handles);
    view3D_calcAxis(1);
    view3D_handles = guidata(view3D_handles.view3DGUI);
end

view3D_handles.arrowHandles = arrowHandles;
view3D_handles.colorList = colorList;
view3D_handles.nspots = nspots;
view3D_handles.goodTime = goodTime;

%save colormap properties so that other programs don't need to reload that
%all the time
view3D_handles.nTags = size(idlist(t1).linklist,1);
view3D_handles.cMap = cMap;
view3D_handles.cMapFact = cMapFact;
view3D_handles.labelList = idlist(1).stats.labelcolor;

%clear lineH
view3D_handles.lineH = [];

%save data
guidata(view3D_handles.view3DGUI,view3D_handles);

close(waitbarHandle);
view3D_refresh;
