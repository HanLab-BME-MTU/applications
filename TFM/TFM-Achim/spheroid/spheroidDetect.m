function [maskSphrd,edgeSphrd,cellDistFromEdge,distGrad,toDoList]=spheroidDetect(fileListSpheroid,fileListMonolayer,smoothing,toDoList,Pix,option,doSave)
clear persistent

try
    load('fileAndFolderNames');
catch
    display('Browse to the project folder!')
    return;
end

%read in Stack of images:
if nargin < 1 || isempty(fileListSpheroid)
    try
        fileListSpheroid =getFileListFromFolder(['data',filesep,'6SpheroidFinal']);
    catch
        [filename, pathname] = uigetfile({'*.tif';'*.jpg';'*.png';'*.*'}, ...
            'Select First Spheroid Image');
        
        if ~ischar(filename) || ~ischar(pathname)
            return;
        end
        
        fileListSpheroid = getFileStackNames([pathname filesep filename]);
    end
else
    isValid = 1;
    for frame = 1:numel(fileListSpheroid)
        isValid = isValid && exist(fileListSpheroid{frame}, 'file');
    end
    if ~isValid
        error('Invalid input files.');
    end
end

if nargin < 2 || isempty(fileListMonolayer)
    try
         fileListMonolayer=getFileListFromFolder(['data',filesep,'6MonolayerFinal']);
    catch
       [filename, pathname] = uigetfile({'*.tif';'*.jpg';'*.png';'*.*'}, ...
           'Select First Spheroid Image');

       if ~ischar(filename) || ~ischar(pathname)
           return;
       end

       fileListMonolayer = getFileStackNames([pathname filesep filename]);
    end
else
    isValid = 1;
    for frame = 1:numel(fileListMonolayer)
        isValid = isValid && exist(fileListMonolayer{frame}, 'file');
    end
    if ~isValid
        error('Invalid input files.');
    end
end


if nargin < 3 || isempty(smoothing)
    smoothing=10;
end

if nargin < 4 || isempty(toDoList)
    toDoList=length(fileListSpheroid);
end

if nargin < 5 || isempty(Pix)
    existPix=0;
else
    existPix=1;
end

% This is used to reset e.g. negative intensity values back to background
% levels:
if nargin < 6 || isempty(option)
    option='noholes'; % The bgVal is set to the intensity value at 1% of the intensity spectrum.
end

if nargin < 7 || isempty(doSave)
    doSave=1; % The bgVal is set to the intensity value at 1% of the intensity spectrum.
end

doPlot=0;

for frame=toDoList
    text=['Detect sheet edges in ',num2str(toDoList(end)),' images'];
    progressText(frame/toDoList(end),text);
    
    currImageSphrd=double(imread(fileListSpheroid{frame}));    
    currImageMonol=double(imread(fileListMonolayer{frame}));    
    [rowsOrg , colsOrg ]=size(currImageSphrd);
     
    %**********************************************************************
    % 1.Step: crop non-zero part of the image                             *
    %**********************************************************************
    if ~existPix
        % Crop an inner part of the image, the boundary might contain zero
        % pixels. The following lines find the largest rectangle with non-zero
        % values that fits into the image (with equal spacing to the bundary).
        realIm=(currImageSphrd~=0);
        bwBD=bwboundaries(realIm,8,'noholes');
        % one should first finde the maximum, anyways:
        bwBD=bwBD{1};   

        if doPlot==1
            figure, imshow(currImageSphrd,[]), title('Croped part of the image')
            colormap gray;
            hold on
            plot(bwBD(:,2),bwBD(:,1),'*b')
            hold off;
        end

        distVals=[bwBD(:,1),bwBD(:,2),rowsOrg-bwBD(:,1),colsOrg-bwBD(:,2)];
        Pix(frame)=max(min(distVals,[],2));
    end

    IcropSphrd =currImageSphrd((1+Pix(frame)):(rowsOrg-Pix(frame)),(1+Pix(frame)):(colsOrg-Pix(frame)));
    IcropMonol =currImageMonol((1+Pix(frame)):(rowsOrg-Pix(frame)),(1+Pix(frame)):(colsOrg-Pix(frame)));
        
    [rows,cols]=size(IcropSphrd);
    
    %**********************************************************************
    % 2. Subtract background                                              *
    %**********************************************************************
    filterNoise=1;
    %bg is intentionally not substracted!
    %remove noise by filtering image with a Gaussian whose sigma = 10 pixel
    if filterNoise
        ISphrd = filterGauss2D(IcropSphrd,smoothing);
        IMonol = filterGauss2D(IcropMonol,smoothing);
    else
        ISphrd = IcropSphrd;
        IMonol = IcropMonol;
    end

    %**********************************************************************
    % 2. Normalize the images and take their difference                   *
    %**********************************************************************
    %get minumum and maximum pixel values in image
    minSignal = min(ISphrd(:));
    maxSignal = max(ISphrd(:));

    %normalize nonzero value between 0 and 1
    ISphrdNorm= (ISphrd - minSignal) / (maxSignal - minSignal);
    
    %get minumum and maximum pixel values in image
    minSignal = min(IMonol(:));
    maxSignal = max(IMonol(:));

    %normalize nonzero value between 0 and 1
    IMonolNorm= (IMonol - minSignal) / (maxSignal - minSignal);
    
    % add the two images:
    I=IMonolNorm-ISphrdNorm;
    
    %**********************************************************************
    % 3. Step: cut off, prepare output                                    *
    %**********************************************************************
    
    %estimate the bg-level
    [~, levelSphrd]=cutFirstHistMode(ISphrdNorm,0);
   %[~, levelMonol]=cutFirstHistMode(IMonolNorm,0);
    
    % In this case, the mask will only be 1 where there is spheroid signal
    % but no monolayersignal:
    cutOffLvl=-levelSphrd;
  
    bwMask=ones(size(I));
    bwMask(I>cutOffLvl)=0;
    
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
    
    % Remove super rough features along the edge:    
    % Maybe distance transform might work better here!
    se   = strel('disk', smoothing);
    bwMask=imopen(bwMask,se);
    if doPlot==1
        figure, imagesc(bwMask), title('Mask smooth')
    end
    
    % get largest (non-zero) connected component:
    % [~,maskLabeled]=bwboundaries(bwMask,4);
    [cc]=bwconncomp(bwMask,4);
    numPixelsCC = cellfun(@numel,cc.PixelIdxList);
    [~,idx] = max(numPixelsCC);    
    bwMaskFinal=zeros(size(bwMask));
    bwMaskFinal(cc.PixelIdxList{idx})=1;
    
    if doPlot==1
        figure, imagesc(bwMaskFinal), title('B/W mask')
    end
    
    %**********************************************************************
    % 4. Post processing:
    %**********************************************************************
    
    % extract the boundary curve:
    [edgeCurve]=bwboundaries(bwMaskFinal,4);

    % look in woundEdgeDetect in case holes in the spheroid region should
    % NOT be considered as still being part of the monolayer.
    if length(edgeCurve)>1
        display('There are holes in the spheroid foot print!')
    end
    
    edgeCurveFinal=edgeCurve{1};
    edgeCurveFinal(edgeCurveFinal(:,1)==1 | edgeCurveFinal(:,2)==1 | edgeCurveFinal(:,1)==rows | edgeCurveFinal(:,2)==cols,:)=[];
    
    
    if doPlot==1
        figure, imagesc(bwMaskFinal), title('Mask smooth')
        hold on
        plot(edgeCurveFinal(:,2) ,edgeCurveFinal(:,1) ,'ko');
        hold off
    end%         figure(2), imagesc(maskSphrd(frame).mat), title('Mask Final')
%         hold on
%         plot(edgeSphrd(frame).pos(:,2) ,edgeSphrd(frame).pos(:,1) ,'ko');
%         hold off

    
    %**********************************************************************
    % 5. Prepare the output:
    %**********************************************************************
    
    % shift to the right coordinate system.
    edgeSphrd(frame).pos = edgeCurveFinal +Pix(frame);
    
    % shift bwMaskFinal
    [rpos,cpos]=ind2sub(size(bwMaskFinal),find(bwMaskFinal==1));
    rpos=rpos+Pix(frame);
    cpos=cpos+Pix(frame);

    maskSphrd(frame).mat=false(size(currImageSphrd));
    linInd=sub2ind(size(maskSphrd(frame).mat),rpos,cpos);
    maskSphrd(frame).mat(linInd)=true;    
    
    % plot the final results:
    if doPlot==1
        figure(1)
        subplot(1,2,1);
        imagesc(currImageMonol), title(['Final Results',num2str(frame)]);
        hold on;
        plot(edgeSphrd(frame).pos(:,2) ,edgeSphrd(frame).pos(:,1) ,'-r');        
        colormap gray;
        hold off
        
        subplot(1,2,2);
        imagesc(currImageSphrd), title('Final Results');
        hold on;
        plot(edgeSphrd(frame).pos(:,2) ,edgeSphrd(frame).pos(:,1) ,'-r');        
        colormap gray;
        hold off

%         figure(2), imagesc(maskSphrd(frame).mat), title('Mask Final')
%         hold on
%         plot(edgeSphrd(frame).pos(:,2) ,edgeSphrd(frame).pos(:,1) ,'ko');
%         hold off
    end
    
    % now bring the curve in x-y coordinates, this is different from
    % woundEdgeDetect!:
    edgeSphrd(frame).pos=fliplr(edgeSphrd(frame).pos);
end

if doSave
    save([path_mechTFM,filesep,'xSpheroidEdge.mat'],'maskSphrd','edgeSphrd');
end

if nargout>3
    r=5;
    [cellDistFromEdge,distGrad]=createCellDistMat(maskSphrd,r,Pix,toDoList);
    if doPlot==1
        figure;
        frame=toDoList(end);
        imagesc(cellDistFromEdge(:,frame).mat), title('cell distance from edge');
        colormap jet;
    end
end