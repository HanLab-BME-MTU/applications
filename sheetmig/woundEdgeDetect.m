function [sheetMask,sheetBnD,sheetEdge,cellDistFromEdge,distGrad,toDoList]=woundEdgeDetect(imageFileList,r,smoothFac,edgeFilter,sigma,Pix,option)
%read in Stack of images:
if nargin < 1 || isempty(imageFileList)
   [filename, pathname] = uigetfile({'*.tif';'*.jpg';'*.png';'*.*'}, ...
       'Select First Phase Contrast Image');
   
   if ~ischar(filename) || ~ischar(pathname)
       return;
   end
   
   imageFileList = getFileStackNames([pathname filesep filename]);
else
    isValid = 1;
    for frame = 1:numel(imageFileList)
        isValid = isValid && exist(imageFileList{frame}, 'file');
    end
    if ~isValid
        error('Invalid input files.');
    end
end

if nargin<2 || isempty(r)
    r=5;
end
% make sure that r is odd:
r=r+(~mod(r,2));

if nargin < 3 || isempty(smoothFac)
    smoothFac=2;
end

if nargin < 4 || isempty(edgeFilter)
    edgeFilter='sobel';
elseif ~strcmp(edgeFilter,'sobel') && ~strcmp(edgeFilter,'canny') && ~strcmp(edgeFilter,'prewitt') && ~strcmp(edgeFilter,'none')
    display('This filter is not supported');
    return;
end

% This is only used for the canny filter:
if nargin < 5 || isempty(sigma)
    sigma=2;
end

if nargin < 6 || isempty(Pix)
    existPix=0;
else
    existPix=1;
end

% This is used to reset e.g. negative intensity values back to background
% levels:
if nargin < 7 || isempty(option)
    option='noholes'; % The bgVal is set to the intensity value at 1% of the intensity spectrum.
end

doPlot=0;

% All image files need to be precessed. If the cell sheet fills up the
% whole field of view (or some other problem occurs), then, this list is
% shortened:
toDoList=1:length(imageFileList);

for frame=toDoList
    text=['Detect sheet edges in ',num2str(toDoList(end)),' images'];
    progressText(frame/toDoList(end),text);
    
    currImage=double(imread(imageFileList{frame}));    
    [rowsOrg, colsOrg]=size(currImage);
    
    %**********************************************************************
    % 1.Step: crop non-zero part of the image                             *
    %**********************************************************************
    if ~existPix
        % Crop an inner part of the image, the boundary might contain zero
        % pixels. The following lines find the largest rectangle with non-zero
        % values that fits into the image (with equal spacing to the bundary).
        realIm=(currImage~=0);
        bwBD=bwboundaries(realIm,8,'noholes');
        % one should first finde the maximum, anyways:
        bwBD=bwBD{1};   

        if doPlot==1
            figure, imshow(currImage,[]), title('Gradient magnitude (gradmag)')
            colormap gray;
            hold on
            plot(bwBD(:,2),bwBD(:,1),'*b')
            hold off;
        end

        distVals=[bwBD(:,1),bwBD(:,2),rowsOrg-bwBD(:,1),colsOrg-bwBD(:,2)];
        Pix(frame)=max(min(distVals,[],2));
    end

    Icrop =currImage((1+Pix(frame)):(rowsOrg-Pix(frame)),(1+Pix(frame)):(colsOrg-Pix(frame)));
    I=Icrop;
    
    [rows,cols]=size(Icrop);
    
    %**********************************************************************
    % 2. Step: Preprocessing: gradient subtraction, doublelog, G-filter   *
    %**********************************************************************
    
    % calculate the gradient information in the image:
    if strcmp(edgeFilter,'sobel') || strcmp(edgeFilter,'prewitt')
        hy = fspecial(edgeFilter); %'prewitt' % perform very much the same!
        hx = hy';
        Iy = imfilter(double(I), hy, 'replicate');
        Ix = imfilter(double(I), hx, 'replicate');
        gradmag = sqrt(Ix.^2 + Iy.^2);
    elseif strcmp(edgeFilter,'canny')
        [gradmag] = steerableFiltering(I,1,sigma);
    elseif strcmp(edgeFilter,'none')
        gradmag=zeros(size(I));
    end
    
    if doPlot==1
        figure, imshow(gradmag,[]), title('Gradient magnitude')
        colormap gray;
    end
    
    se   = strel('disk', r);
    Gc = imclose(gradmag,se);
    if doPlot==1
        figure, imshow(Gc,[]), title('Gradient image closed')
        colormap gray;
    end
        
    % cut off the BG:
    [~, level]=cutFirstHistMode(Gc,0);
    bwMask=Gc>level;
    
    if doPlot==1
        figure, imagesc(bwMask), title('B/W mask')
        colormap gray;
    end
    
    % Fill up the holes in the regions if wanted:
    if strcmp(option,'noholes');
        bwMask=imfill(bwMask,'holes');
        if doPlot==1;
            figure, imagesc(bwMask), title('B/W mask')
        end
    end
     
    % get largest (non-zero) connected component:
    [~,maskLabeled]=bwboundaries(bwMask,4);
        
    if doPlot==1
        figure, imagesc(maskLabeled), title('B/W mask')
    end
    
    histL=maskLabeled(bwMask>0);
    n = histc(histL(:),min(histL(:)):max(histL(:)));
    % find the two larges regions (either the sheet or the bg):
    LofMax=find(n==max(n));
    
    maskSheet=false(size(maskLabeled));
    maskSheet(maskLabeled==LofMax)=true;
    
    if doPlot==1
        figure, imagesc(maskSheet), title('Mask for the cell sheet')
    end
    
    coveredAreaRatio=sum(maskSheet(:))/numel(maskSheet);
    if coveredAreaRatio>0.90 || coveredAreaRatio<0.1
        display(['Max area ratio exceeded. Most likely cells filled up the whole field of view at frame: ',num2str(frame)]);
        toDoList=1:frame-1;
        break;
    end
    % Post processing:
    % Remove super rough features along the edge:    
    % Maybe distance transform might work better here!
    se   = strel('disk', smoothFac*r);
    maskSmooth=imopen(maskSheet,se);
    if doPlot==1
        figure, imagesc(maskSmooth), title('Mask smooth')
    end
    
    % extract the boundary curve:
    [~,maskSmoothLbl]=bwboundaries(maskSmooth,4);

    histSmoothL=maskSmoothLbl(maskSmooth>0);
    n = histc(histSmoothL(:),min(histSmoothL(:)):max(histSmoothL(:)));
    % find the largest region, the sheet:
    LofSmoothMax=find(n==max(n));    
    % But this still might contain holes at the boundary:
    maskWithBndHoles=false(size(maskSmoothLbl));
    maskWithBndHoles(maskSmoothLbl==LofSmoothMax)=true;
    
    % The largest hole is the background:
    [~,maskHolesLbl]=bwboundaries(~maskWithBndHoles,4);
    listHolesLbl=maskHolesLbl(~maskWithBndHoles>0);
    nHoles = histc(listHolesLbl(:),min(listHolesLbl(:)):max(listHolesLbl(:)));
    % find the holes, (which are smaller then bg, the largest hole):
    [LofHoles]=find(nHoles~=max(nHoles));
    
    maskFinal=maskWithBndHoles;
    % Now fill up the wholes:
    for k=1:length(LofHoles)
        maskFinal(maskHolesLbl==LofHoles(k))=true;
    end
    
    % now get again the final boundaries:
    [bndCurvesFinalStr]=bwboundaries(maskFinal,4);
    % There should be only one boundary left:
    if length(bndCurvesFinalStr)==1
        bndCurveFinal=bndCurvesFinalStr{1};
    else
        error('Something went wrong there is more than one domain left!')
    end
    
    edgeCurveFinal=bndCurveFinal;
    edgeCurveFinal(edgeCurveFinal(:,1)==1 | edgeCurveFinal(:,2)==1 | edgeCurveFinal(:,1)==rows | edgeCurveFinal(:,2)==cols,:)=[];
    
    
    if doPlot==1
        figure, imagesc(maskFinal), title('Mask smooth')
        hold on
        plot(bndCurveFinal(:,2) ,bndCurveFinal(:,1) ,'b-');
        plot(edgeCurveFinal(:,2),edgeCurveFinal(:,1),'ro');
        hold off
    end
    
    % check if the bndCurveFinal touches all boundaries. Then stop
    % analysis:
    [rmax,cmax]=size(maskFinal);
    if max(bndCurveFinal(:,1))==rmax && min(bndCurveFinal(:,1))==1 && max(bndCurveFinal(:,2))==cmax && min(bndCurveFinal(:,2))==1
        display(['The cell sheet reached the field of view boundaries in frame: ',num2str(frame)]);
        toDoList=1:frame-1;
        break;
    end
    
    % prepare the output:
    
    % shift to the right cordinate system.
    sheetBnD(frame).pos = bndCurveFinal +Pix(frame);
    sheetEdge(frame).pos= edgeCurveFinal+Pix(frame);
    
    % Sometimes it happens that the first and last points are the same:
    %if compPts(sheetEdge(frame).pos(1,:),sheetEdge(frame).pos(end,:))
    %    if sheetEdge(frame).pos(1,1)<rmax/2
    %        sheetEdge(frame).pos(end,:)=[];
    %    else
    %        sheetEdge(frame).pos(1,:)=[];
    %    end
    %end

    
    % shift maskFinal
    [rpos,cpos]=ind2sub(size(maskFinal),find(maskFinal==1));
    rpos=rpos+Pix(frame);
    cpos=cpos+Pix(frame);

    sheetMask(frame).mat=false(size(currImage));
    linInd=sub2ind(size(sheetMask(frame).mat),rpos,cpos);
    sheetMask(frame).mat(linInd)=true;    
    
    % plot the final results:
    if doPlot==1
        figure;
        Idspl=currImage;
        Idspl(sheetMask(frame).mat == 1) = 0;
        imagesc(Idspl), title('Extended nuclei');
        colormap jet;
    end

    % plot the final results:
%     figure(1)
%     imagesc(currImage), title('Extended nuclei');
%     hold on
%         plot(sheetBnD(frame).pos(:,2) ,sheetBnD(frame).pos(:,1) ,'k-');
%         plot(sheetEdge(frame).pos(:,2),sheetEdge(frame).pos(:,1),'r.');
%     hold off
%     colormap jet;
end

if nargout>3    
    [cellDistFromEdge,distGrad]=createCellDistMat(sheetMask,r,Pix,toDoList);
    if doPlot==1
        figure;
        frame=toDoList(end);
        imagesc(cellDistFromEdge(:,frame).mat), title('cell distance from edge');
        colormap jet;
    end
end