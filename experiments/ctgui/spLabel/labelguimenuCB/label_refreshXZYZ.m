function label_refreshXZYZ
% function to update xz/yz maximum projection of 3D movie stacks in labelgui

%get necessary handles and data

%zFigureData
imgFigH = GetUserData(openfig('labelgui','reuse'),'currentWindow');
zFigH = GetUserData(imgFigH,'XZYZFigureH');

%labelguiH
labelguiH = findall(0, 'Tag', 'labelgui');
labelguiHandles = guidata(labelguiH);

if isempty(zFigH) %some user has closed figure window. Reopen
    hObject=labelguiHandles.label_showXZYZ;
    zFigH = label_showXZYZ_CB(hObject,[],[],1);
end

axesH = GetUserData(zFigH, 'axesH');
blowUpVector = GetUserData(zFigH, 'blowUpVector');
pixelRatio = GetUserData(zFigH, 'pixelRatio');
pixelSize = GetUserData(zFigH, 'pixelSize');
movieSize = GetUserData(zFigH, 'movieSize');

%labelPanelData
movie = GetUserData(imgFigH,'mv'); 
idlist = GetUserData(imgFigH,'idlist');
dataProperties = GetUserData(imgFigH,'dataProperties');
PIXELSIZE_XY=dataProperties.PIXELSIZE_XY;
PIXELSIZE_Z=dataProperties.PIXELSIZE_Z;

%labelguiData
curr_time = round(get(labelguiHandles.slider3,'Value'));
if ~isempty(idlist)
    cMap = GetUserData(labelguiH,'cMap');
    cMapFact = size(cMap,1)/idlist(1).stats.maxColor;
end

sliceMode = get(labelguiHandles.radiobutton2, 'Value');
if sliceMode
    curr_slice = get(labelguiHandles.slider4, 'Value');
else
    curr_slice = [];
end

bigPicture = movie(:,:,blowUpVector,1,curr_time);

%check wheter we use true brightness or whether the max intensity of
    %this image is 1
    useTrueBrightness = strcmp(get(labelguiHandles.label_trueBrightness,'Checked'),'on');
    if ~useTrueBrightness
        bigPicture = bigPicture/max(bigPicture(:));
    end

%show data
for xy=1:2
    
    %get projections
    projection = max(bigPicture,[],xy); %the horizontal axis in the image is the first dimesion in the movie!
    
    %make z the rows and x or y the cols and remove third (singleton) dimension
    projection = squeeze(projection)';
    
    if xy==2 %turn img around so that first view is from bottom and second from the right
        projection = projection(:,end:-1:1);
    end
    
    
    %show projections
    axes(axesH(xy));
    cla;
    imshow(projection);
    hold on;
    
    if ~isempty(idlist)
        %---------------------plot labels------------------------------------------
        
        if ~isempty(idlist(curr_time).linklist)
            for i=1:max(idlist(curr_time).linklist(:,2))
                
                %find colors belonging to same spot
                rowIdx=find(idlist(curr_time).linklist(:,2)==i);
                
                %read and change coordinates
                coord=idlist(curr_time).linklist(rowIdx(1),9:11)./[PIXELSIZE_XY PIXELSIZE_XY PIXELSIZE_Z];
                coord(3) = (movieSize(3) - coord(3)+0.75)*pixelRatio;
                coord(2) = movieSize(2)-coord(2)+1;
                
                %list of colors in one spot
                cList=idlist(curr_time).linklist(rowIdx,4);
                
                
                if length(cList)==1
                    % plot position (xy-z)
                    pH=plot(coord(xy),coord(3),'o','MarkerSize',4);
                    labelColor(1,:)=cMap(cList*cMapFact,:);
                    set(pH,'MarkerEdgeColor',labelColor);
                    type=char(idlist(1).stats.labelcolor(log2(cList)+1));
                else
                    %plot all spots making up this one
                    for j=length(cList):-1:1
                        pH=plot(coord(xy),coord(3),'o','MarkerSize',4+4*(j-1));
                        labelColor(j,:)=cMap(cList(j)*cMapFact,:);
                        set(pH,'MarkerEdgeColor',labelColor(j,:));
                    end
                    labelColor=mean(labelColor,1);
                    type='fus';
                end
                
                % write label
                txtH=text(coord(xy)+5,coord(3),type,'Color',labelColor);
                %set label2spot number
                set(txtH,'UserData',i);
                % popup CB
                set(txtH,'ButtonDownFcn','popuphere');    
            end;
            
        end;
        %---------------------------------------------------------------------------
    end
    
    
    %plot current slice (draw rectangle)
    if sliceMode
        %draw rectangle (x y w h)
        rH = rectangle('Position',[0.5,(movieSize(3)-curr_slice)*pixelRatio+0.5,movieSize(xy),pixelRatio],'EdgeColor','r');
    end
    
end %for-loop

%make both images the same size
axis(axesH,'normal');

% make gui top most
figure(imgFigH);
figure(labelguiH);