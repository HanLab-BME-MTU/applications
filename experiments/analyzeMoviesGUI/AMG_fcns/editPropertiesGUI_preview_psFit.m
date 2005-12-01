function editPropertiesGUI_preview_psFit(hObject,eventdata,handles)
%Callback for button appearing in the editPropertiesGUI-previewWindow for the
%PSF-fitting

toggleState = get(hObject,'Value');

if toggleState == 0
    %button is being untoggled
    
    %could mean cancel all in the future, but not yet
    set(hObject,'Value',1);
    return
end

%-------------check that there are any valid points

%read necessary data
previewFigH = handles.previewPSF.figureH;
coords = handles.previewPSF.cord;
frameList = handles.previewPSF.selFrames;
axesList = handles.previewPSF.axesList;
imgSize = handles.previewPSF.imgSize;
%spots have to be at least 7 sigma apart
minDist = ceil(7*handles.FILTERPRM(1));
%spots have to be at least 3.5 sigma from the border
minBorderDist = ceil(3.5*handles.FILTERPRM(1:3));
%use 1.5*sigma, rounded up, for patchSize
deltaXYZ = ceil([handles.FILTERPRM(1),handles.FILTERPRM(2),handles.FILTERPRM(3)]*1.5);
oldSigma = handles.FILTERPRM([1,3]);
oldSigmaCorr = handles.sigmaCorrection([1,2]);

for i = 1:length(axesList)
    %for all images: if there are points that are at least 2*3*sigmaXY pixels distant from
    %each other, make them visible; delete other entries in coordList
    %since we can only cut squares and not circles, the minDist criterion has to
    %be met either one or the dimension
    
    coordList(i).coord = cat(1,coords(i).sp.cord);
    nPoints = size(coordList(i).coord,1);
    
    dmX = distMat(coordList(i).coord(:,1));
    dmY = distMat(coordList(i).coord(:,2));
    dia = diag(ones(nPoints,1)*999);
    
    %there are only zeros in the diagonal -> correct!
    minX = min(dmX+dia,[],2);
    minY = min(dmY+dia,[],2);
    
    %remove entries form coordList that lie too close to neighbouring points
    %or too close to any of the image edges
    rmIdx = find( (minX<minDist & minY<minDist) |...
        sum(floor(coordList(i).coord-repmat(minBorderDist,nPoints,1)) < repmat([0 0 0],nPoints,1),2) | ...
        sum(ceil(coordList(i).coord+repmat(minBorderDist,nPoints,1)) > repmat(imgSize([2,1,3]),nPoints,1),2));
    
    
    coordList(i).coord(rmIdx,:) = [];
    
    %plot green circles around good spots
    axes(axesList(i));
    if ~isempty(coordList(i).coord)
        plot(coordList(i).coord(:,1),coordList(i).coord(:,2),'og');
    end
end

%make sure there are any spots left
if isempty(cat(1,coordList.coord))
    h = errordlg('Sorry, there are no valid signals to fit the psf (too close together)','Try again');
    uiwait(h)
    set(hObject, 'Value', 0);    
    return
end

%--------------------------------------------------

%tell the user what he/she's doing
ans = questdlg(['You are about to start fitting the psf to the microscopy data. To this end', ...
    'please draw rectangles around 2-3 isolated tag signals. If you plan to do a background subtraction and haven''t done so yet',...
    ' cancel now and apply the respective settings first.'],'Please follow the instructions carefully','Continue','Cancel','Continue');

if ~strcmp(ans,'Continue')
    %user cancelled
    set(hObject, 'Value', 0);
    return
end



enough = 0;
selectedRectangleH = [];
psFitData = struct('patchXYZ',[]);

while ~enough
    
%     if ~isempty(selectedRectangleH)
%         delete(selectedRectangleH)
%     end
    
    h = helpdlg('Please select a signal (click on image first)');
    uiwait(h)
    
    %this part is from help rbbox: draw rectangle
    k = waitforbuttonpress; %make axes current
    k = waitforbuttonpress;
    point1 = get(gca,'CurrentPoint');    % button down detected
    finalRect = rbbox;                   % return figure units
    point2 = get(gca,'CurrentPoint');    % button up detected
    point1 = point1(1,1:2);              % extract x and y
    point2 = point2(1,1:2);
    selectPosMin = min(point1,point2);        % lower left corner
    selectPosMax = max(point1,point2);        % upper right corner
    
    if ~ishandle(previewFigH) %user cancelled by closing previewWindow
        delete(gca); %calling gca above created axes on the EPGUI-window
        return %end evaluation here
    end
    
    if any(selectPosMin == selectPosMax) | any([selectPosMin,selectPosMax]<0)
        %user selecte only a point/line, or has not selected an image (in that case, the lower left most axes are current)
        h = warndlg(['Please click first on the image and then select a good signal', ...
            '(green circle) by drawing a recangle around it'],'Insufficient input');
        uiwait(h);
        
    else %good input
        
        %find the window in which the selection has been made
        currentAxes = gca;
        currentIdx = find(axesList==currentAxes);
        
        currentCoords = coordList(currentIdx).coord;
        nCoords = size(currentCoords,1);
        
        %find the selected points that are within the rectangle, i.e. lie between selectPosMin and selectPosMax
        goodIdx = find(sum([currentCoords(:,1:2)>repmat(selectPosMin,nCoords,1),currentCoords(:,1:2)<repmat(selectPosMax,nCoords,1)],2)==4);
        
        if ~isempty(goodIdx) & length(goodIdx)<2
            %calculate patchSize (has to be odd size!)
            %if center is whithin pixel 4, we want to take pixels 1:7 -> floor -> +-same distance
            patchXYZ = [floor(currentCoords(goodIdx,:)-deltaXYZ);floor(currentCoords(goodIdx,:)+deltaXYZ)];
            rectPosMin = patchXYZ(1,1:2)+[-0.5 0.5];
            rectPosDelta = diff(patchXYZ(:,1:2),1,1)+1;
            
            %change patchXYZ into matrix coordinates and store
            psFitData(end+1).patchXYZ = patchXYZ(:,[2,1,3]);
            psFitData(end).imgNumber = frameList(currentIdx);
            psFitData(end).relSignalPosition = currentCoords(goodIdx,[2,1,3])-psFitData(end).patchXYZ(1,:);
           
            
            %draw rectangle. make it the correct size. correct for rectangle properties
            selectedRectangleH = rectangle('Position',[rectPosMin,rectPosDelta],'EdgeColor','r');
            figure(previewFigH);
            
            again = questdlg('Do you want to select another signal?','','Yes','No','Cancel','Yes');
            if ~strcmp(again,'Yes')
                enough = 1;
            end
            if strcmp(again,'Cancel')
                set(hObject, 'Value', 0);
                return
            end
        else
            %either two spots or none
            h = helpdlg(['You have selected ',num2str(length(goodIdx)),' signals. Try again'],'Not quite good enough');
            uiwait(h);
        end %if ~isempty(goodIdx) & length(goodIdx)<2
        
    end %test if rectangle. if any(selectPosMin == selectPosMax)
    
end % while ~enough

%make sure there is at least one signal selected
if isempty(psFitData)
    errordlg('no signal selected');
    set(hObject, 'Value', 0);
    return
end

psFitData(1) = [];
imgNumberList = cat(1,psFitData.imgNumber);
nImg = length(imgNumberList);

%load all images
switch handles.previewData.movieName(end)
    case 'c'
        imgList = readmat(handles.previewData.movieName,imgNumberList);
    case 'r'
        imgList = zeros(handles.header.numRows, handles.header.numCols, handles.header.numZSlices, 1, nImg);
        for i = 1:nImg
            imgList(:,:,:,:,i) = r3dread(handles.previewData.movieName,imgNumberList(i),1);
        end
end

%do fitting
for i = 1:nImg 
    
    %get current Img (for bg-subtraction)
    currentImg = squeeze(imgList(:,:,:,:,i));
    
    %get psf-Img from raw data
    psfImg = currentImg(psFitData(i).patchXYZ(1,1):psFitData(i).patchXYZ(2,1),psFitData(i).patchXYZ(1,2):psFitData(i).patchXYZ(2,2),...
        psFitData(i).patchXYZ(1,3):psFitData(i).patchXYZ(2,3));
    
    %test noise
    for j = 1:imgSize(3)
        nse(j) = imNoiseEstim(currentImg(:,:,j));
    end
    currentNoise = mean(nse);
    
    spsf = (size(psfImg)+1)/2;
    %lower bound for SNR
    currentMinSNR = psfImg(spsf(1),spsf(2),spsf(3))/currentNoise-2
    
%     %background: 2-pixels area around the patch
%     bgImg1 = currentImg([psFitData(1).patchXYZ(1,1)-2:psFitData(1).patchXYZ(1,1)-1,psFitData(1).patchXYZ(2,1)+1:psFitData(1).patchXYZ(2,1)+2],...
%         psFitData(i).patchXYZ(1,2):psFitData(i).patchXYZ(2,2), psFitData(i).patchXYZ(1,3):psFitData(i).patchXYZ(2,3));
%     bgImg2 = currentImg(psFitData(i).patchXYZ(1,1):psFitData(i).patchXYZ(2,1),...
%         [psFitData(1).patchXYZ(1,2)-2:psFitData(1).patchXYZ(1,2)-1,psFitData(1).patchXYZ(2,2)+1,psFitData(1).patchXYZ(2,2)+2],...
%         psFitData(i).patchXYZ(1,3):psFitData(i).patchXYZ(2,3));
%     bgImg3 = currentImg(psFitData(i).patchXYZ(1,1):psFitData(i).patchXYZ(2,1),psFitData(i).patchXYZ(1,2):psFitData(i).patchXYZ(2,2),...
%         [psFitData(1).patchXYZ(1,3)-2:psFitData(1).patchXYZ(1,3)-1,psFitData(1).patchXYZ(2,3)+1:psFitData(1).patchXYZ(2,3)+2]);
%     bg = mean([bgImg1(:);bgImg2(:);bgImg3(:)]);
    
 %    psm = mean(psfImg(:));
 imm = mean(currentImg(:))
 imm2 = median(currentImg(:))
%  imm2 = median(currentImg(:));
%  max([psm,imm]);
%  psfImg1 = psfImg-max([psm,imm]);
% %     
%     psfImg2 = psfImg - bg*1.2;
% %     psfImg3 = psfImg - imm2;
% psfImg3 = psfImg-bg;
    
    %prepare initial parameters (dx, dy, dz, amp (not used in the fit), sigXY, sigZ)
    parms1 = [psFitData(i).relSignalPosition-size(psfImg)/2, 0, oldSigma(1), oldSigma(2)];
    
    [modelPsf,fittedParms1(i,:), mA]=fitGauss2PSF(psfImg,parms1,imm);
    [modelPsf,fittedParms2(i,:), mA]=fitGauss2PSF(psfImg,parms1,imm2);
end

fittedParms1
fittedParms2
oldSigmaCorr
sigmaCorrection2 = median(fittedParms2(:,[5,6]),1).*oldSigmaCorr./oldSigma