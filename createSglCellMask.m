function [sheetMask,sheetBnD,sheetEdge,cellDistFromEdge,distGrad,toDoList]=createSglCellMask(imageFileList,r,Pix)
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

if nargin < 2 || isempty(r)
    r=1;
end

if nargin < 3 || isempty(Pix)
    existPix=0;
else
    existPix=1;
end

% All image files need to be precessed. If the cell sheet fills up the
% whole field of view (or some other problem occurs), then, this list is
% shortened:
toDoList=1:length(imageFileList);

doPlot=1;

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

    sheetMask(frame).mat=false(size(currImage));
    sheetMask(frame).mat((1+Pix(frame)):(rowsOrg-Pix(frame)),(1+Pix(frame)):(colsOrg-Pix(frame)))=true;
    
    % now get the boundaries:
    [bndCurvesStr]=bwboundaries(sheetMask(frame).mat,4);
    % There should be only one boundary left:
    if length(bndCurvesStr)==1
        bndCurve=bndCurvesStr{1};
    else
        error('Something went wrong there is more than one domain left!')
    end
    
    sheetBnD(frame).pos  = bndCurve;
    sheetEdge(frame).pos = [Inf Inf];
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