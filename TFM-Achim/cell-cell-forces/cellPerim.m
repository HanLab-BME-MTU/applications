function [cellEdge]=cellPerim(imageFileList,dilationR,sigmaGauss,closureRadius,holes,thrHoleLen,doPlot,pauseSec)

if nargin <1 || isempty(imageFileList)
   [filename, pathname] = uigetfile({'*.TIF;*.tif;*.jpg;*.png;*.*'}, ...
       'Select first image to be analyzed');
   
   if ~ischar(filename) || ~ischar(pathname)
       return;
   end
   
   imageFileList = getFileStackNames([pathname filesep filename]);
elseif ischar(imageFileList) && ~isdir(imageFileList)
    %then, there is only one file
    dummy=imageFileList;
    clear imageFileList;
    imageFileList{1}=dummy;
else
    isValid = 1;
    for i = 1:numel(imageFileList)
        isValid = isValid && exist(imageFileList{i}, 'file');
    end
    if ~isValid
        error('Invalid input files.');
    end
end

if nargin < 2 || isempty(dilationR)
    dilationR=0;
end

if nargin < 3 || isempty(sigmaGauss)
    sigmaGauss=5; % 3 is the default
end

if nargin < 4 || isempty(closureRadius)
	closureRadius = 3; % 3 is the default
end

if nargin < 5 || isempty(holes) || holes~=1
    holes=0;
end

if (nargin < 6 || isempty(thrHoleLen)) && holes==1
    thrHoleLen=100;
elseif holes==0
    thrHoleLen=-1;
end

if nargin < 7 || isempty(doPlot)
	doPlot = 0;
end

if nargin < 8 || isempty(pauseSec)
	pauseSec = 0;
end


for frameIndex=1:length(imageFileList)
    %read in the image:
    currentImage = double(imread(imageFileList{frameIndex}));    

    % Due to the stage drift correction, the image might be padded with
    % zeros. We have to skip these values. In addition, the first finite
    % values might be bad due to the sub-pixel interpolation. Therefore we
    % introduce a saftyshift of one pixel:
    shift=1;
    
    [rows,cols]=size(currentImage);
    midCol=round(cols/2);
    midRow=round(rows/2);
    
    % with indices I(i,j)
    first_row = find(currentImage(:,midCol),1, 'first')+shift;
    last_row  = find(currentImage(:,midCol),1, 'last')-shift;
    
    first_col = find(currentImage(midRow,:),1, 'first')+shift;
    last_col  = find(currentImage(midRow,:),1, 'last')-shift;

    % The clean image to be analyzed:
    currentImageNonZero=currentImage(first_row:last_row,first_col:last_col);
    
    % Perform the edge detection based on finding a good threshold:
    cellEdge{frameIndex}.mask =zeros(size(currentImage));
    cellEdge{frameIndex}.mask(first_row:last_row,first_col:last_col) =firstMinAfterFirstMaxSeg(currentImageNonZero, closureRadius, sigmaGauss, holes);
    
    % Read out the boundary(ies):
    if holes==1
        B=bwboundaries(cellEdge{frameIndex}.mask,'holes');
    else
        B=bwboundaries(cellEdge{frameIndex}.mask,'noholes');
        % Check if everything is OK. If B>1 then more than one boundary has
        % been detected, e.g. also within the cells. This should not happen.
        if length(B)>1
            display(['There is a problem with image: ',num2str(frameIndex)])
            display(['The number of boundaries is:   ',num2str(length(B))])
        end
    end
    
    % Pick only significant holes based on a threshold value for the
    % length of the hole perimeter:
    for seg=1:length(B)
        segLen(seg)=length(B{seg});
    end
    [segLen,id]=sort(segLen,'descend');
    sigId=id(segLen>thrHoleLen);
    k=1;
    for id=sigId
        Bsig{k}=B{id};
        k=k+1;
    end

    % This should be the perimeter of the whole cell cluster:
    cellEdge{frameIndex}.curve(:,1)=B{1}(:,2);
    cellEdge{frameIndex}.curve(:,2)=B{1}(:,1);
    
    % The remaining entries in the list should be the holes:
    if length(Bsig)>1
        for holeId=2:length(Bsig)
            cellEdge{frameIndex}.hole{holeId-1}.curve(:,1)=Bsig{holeId}(:,2);
            cellEdge{frameIndex}.hole{holeId-1}.curve(:,2)=Bsig{holeId}(:,1);
        end
    end    
   
    % Define the structural element, (the higher dilationR, the smoother 
    % the boundary):
    se = strel('ball', dilationR, 0, 8);
    
    % Dilate the edge a little bit:
    cellEdge{frameIndex}.maskDilated = imdilate(cellEdge{frameIndex}.mask, se);
    
    % read out the new (smoother) boundary(ies):
    % Here connectivity has to be 4, that is edge to edge connection, 
    % whereas 8 would allow also vertex to vertex connections! In the 
    % latter case crossection can be missed in intersecMaskPolygon.
    B=bwboundaries(cellEdge{frameIndex}.maskDilated,4,'noholes');
    % Check if everything is OK. If B>1 then more than one boundary has 
    % been detected, e.g. also within the cells. This should not happen.
    if length(B)>1 
        display(['There is a problem with image: ',num2str(frameIndex)])
        display(['The number of boundaries is:   ',num2str(length(B))])
    end
    cellEdge{frameIndex}.curveDilated(:,1)=B{1}(:,2);
    cellEdge{frameIndex}.curveDilated(:,2)=B{1}(:,1);
    
    
    % create a 0/1 perimeter mask for the dilated version:
    % cellEdge{frameIndex}.maskPerimDilated=bwperim(cellEdge{frameIndex}.maskDilated)
    
    if doPlot==1
        figure(2)
        imagesc(currentImage)
        colormap('gray')
        hold on
        plot(cellEdge{frameIndex}.curve(:,1),cellEdge{frameIndex}.curve(:,2),'b')
        plot(cellEdge{frameIndex}.curveDilated(:,1),cellEdge{frameIndex}.curveDilated(:,2),'g')
        set(gca,'Ydir','reverse')
        title(['Image frame: ',num2str(frameIndex)])
        hold off

        pause(pauseSec)
    end

    cellEdge{frameIndex}.params.dilationR=dilationR;
    cellEdge{frameIndex}.params.sigmaGauss=sigmaGauss;
    cellEdge{frameIndex}.params.closureRadius=closureRadius;
    cellEdge{frameIndex}.params.holes=holes;
    cellEdge{frameIndex}.params.thrHoleLen=thrHoleLen;    
end






