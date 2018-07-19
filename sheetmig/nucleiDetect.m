function [nuclei, movieInfo, dPix, bwBnDMaskOrg, bwLabelsOrg]=nucleiDetect(imageFileList,r,useDblLog,edgeFilter,sigma,p)
%read in Stack of images:
if nargin < 1 || isempty(imageFileList)
   [filename, pathname] = uigetfile({'*.tif';'*.jpg';'*.png';'*.*'}, ...
       'Select First nuclei Image');
   
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

if nargin < 3 || isempty(useDblLog)
    useDblLog=1;
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

% This is used to reset e.g. negative intensity values back to background
% levels:
if nargin < 6 || isempty(p)
    p=0.01; % The bgVal is set to the intensity value at 1% of the intensity spectrum.
end

doPlot=0;

numFrames=length(imageFileList);
for frame=1:numFrames
    text=['Detect nuclei in ',num2str(numFrames),' images'];
    progressText(frame/numFrames,text);
    
    currImage=double(imread(imageFileList{frame}));    
    [rowsOrg, colsOrg]=size(currImage);
    
    %**********************************************************************
    % 1.Step: crop non-zero part of the image                             *
    %**********************************************************************
    
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
    dPix(frame)=max(min(distVals,[],2));

    Icrop =currImage((1+dPix(frame)):(rowsOrg-dPix(frame)),(1+dPix(frame)):(colsOrg-dPix(frame)));
    I=Icrop;
    
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
     
    % Normalize the gradient according to the magnitude of the image intensity:
    gradmag = gradmag*max(I(:))/max(vertcat(eps,gradmag(:)));
    if doPlot==1
        figure, imshow(gradmag,[]), title('Gradient magnitude (gradmag)')
        colormap gray;
    end
       
    minI=min(I(:));
    maxI=max(I(:));
    
    numPix(frame)=numel(I);
    binEdges=linspace(floor(minI),ceil(maxI),1000);
    n = histc(I(:),binEdges);
    
    ncum=cumsum(n);
    binID=find(ncum>p*numPix(frame),1,'first');
    bgVal=binEdges(binID);
    
    % substract the gradient from the original image to enhance the edges:
    ImGrad=I-gradmag;    
    ImGrad(ImGrad<bgVal)=bgVal;  % bring all neg. or very small values back to bg-level!
    
    if doPlot==1
        figure, imagesc(ImGrad), title('Gradient substracted image')
    end
    
    % Equalize neighboring maxima by applying a double-log:
    bgValScaled=log(log(bgVal));
    if useDblLog && isreal(bgValScaled) && bgValScaled>-Inf
        ImGrad=log(log(ImGrad));
    else
        display(['Frame ',num2str(frame),': Couldnt use log(log(.)) to suppress large next to small maxima!']);
    end
    
    % filter with a gaussian:
    Iflt = filterGauss2D(ImGrad,r);
    if doPlot==1
        figure, imagesc(Iflt), title('Gradient substracted and Gauss-filtered image')
        colormap gray;
    end
    
    %**********************************************************************
    % 4. Step: Find the local maxima in the preprocessed image            *
    %**********************************************************************
    
    % Find the local maxima:
    se   = strel('disk', r);
    Imax=locmax2d(Iflt,getnhood(se),1);
    if doPlot==1
        figure, imagesc(Imax), title('Maxima in the image')
        colormap gray;
        figure;
        hist(Imax(Imax(:)>0),1000)
    end
    
    % cut off the maxima in the noise:
    try
        % This only works for dense sheets!
         level1=thresholdFluorescenceImageFewBg(Imax(Imax(:)>0),doPlot);
        % This only works for significant amount of Bg, but is then more reliable than the above method!
        % level1=thresholdFluorescenceImage(Imax(Imax(:)>0),doPlot,1);
    catch
        % This shouldn't be used anymore. Instead, one should use the
        % algorithm above!
        display('!!!switched to cutFirstHistMode!!!')
        [~, level1]=cutFirstHistMode(Imax(Imax(:)>0),0);
    end
    Imax(Imax(:)<level1)=0;
    Imax(Imax(:)>0)=1;
    
    % Show the first set of maxima:
    se   = strel('disk', round(r/4));
    ImaxDil = imdilate(Imax,se);
    if doPlot==1
        Idspl=Icrop;
        Idspl(ImaxDil == 1) = 0;
        figure, imagesc(Idspl), title('Extended detected nuclei')
    end
    
    %**********************************************************************
    % 4. Step: Find a second set of max using a simple Gauss-filter       *
    %**********************************************************************
    
    % Simply filter the original image with a Gaussian, some of these might
    % have been cut-off if they are small in size and have a huge gradient!
    IfltSimple = filterGauss2D(Icrop,r);
    ImaxSimple=locmax2d(IfltSimple,getnhood(se),1);
    [~, level2]=cutFirstHistMode(ImaxSimple(ImaxSimple(:)>0),0);
    ImaxSimple(ImaxSimple(:)<level2)=0;
    ImaxSimple(ImaxSimple(:)>0)=1;
    
    se   = strel('disk', r);
    ImaxComb = imdilate(Imax+ImaxSimple,se);    
    labelsMaxROI = bwlabel(ImaxComb>0, 4);
    % label each Imax. Find doubled values. These maxima have been falsely
    % merged by regiongrowing:
    maxID=find(Imax(:)==1);
    maxLabels=labelsMaxROI(maxID);
    
    n = histc(maxLabels,1:max(maxLabels));
    dblLables=find(n>1);
    falseROI=ismember(labelsMaxROI,dblLables);
    if doPlot==1
        falseROI(ImaxDil==1)=0;
        figure, imagesc(falseROI), title('Falsely fused regions!')
    end
    % Clean Labels:
    labelsMaxROI(falseROI)=0;
    
    % Set the max that we want to keep into the Labels matrix:
    dblLabelInd=ismember(maxLabels,dblLables);
    maxIdKeep=maxID(dblLabelInd);
    
    maxL=max(labelsMaxROI(:));
    % first fill up the Labels and then append at the end:
    newLabelList=[dblLables',(maxL+1):(maxL+length(maxIdKeep)-length(dblLables))];
    for k=1:length(maxIdKeep)
        labelsMaxROI(maxIdKeep(k))=newLabelList(k);
    end
    
    %figure, imagesc(labelsMaxROI), title('Labels')
    
    regions = regionprops(labelsMaxROI, 'centroid');
    
    nucPos=round(vertcat((regions(:).Centroid)));
    ImaxCombClean=zeros(size(Icrop));
    ImaxCombClean(sub2ind(size(Icrop),nucPos(:,2),nucPos(:,1)))=1;
    
    if doPlot==1
        Idspl=Icrop;
        se   = strel('disk', round(r/4));
        ImaxCombCleanDil = imdilate(ImaxCombClean,se);
        Idspl(ImaxCombCleanDil > 0) = 0;
        figure, imagesc(Idspl), title('Extended nuclei')
    end
    
    %**********************************************************************
    % 5. Step: Watershed to get boundaries for the nuclei                 *
    %**********************************************************************
    
    % Use the detected maxima as forground markers of a watershed:
    if nargout>3
        % Find good background markers, this might be difficult in a very
        % densely populated area. Maybe find local intensity minimas.
        
        % I is log(log(I)) or G-filtered image:
        I=Iflt;
        bgm = I<max(level1);
        %The boundary detection for the nuclei can be improved by first
        %finding the BGmin before doing a watershed!
        
        % Find the regional maxima in the gradient image:
        gradmag2 = imimposemin(gradmag, bgm | ImaxCombClean);
        if doPlot==1
            figure, imshow(gradmag2), title('Thresholded opening-closing by reconstruction (bw)')
        end
        % the regions that grow out of the bgm will contain no intensity
        % maximum and will be sorted out later on.

        % compute the watershed:
        Lws = watershed(gradmag2);
        
        % Visualization:
        if doPlot==1
            se   = strel('disk', round(r/4));
            I4 = (I-min(I(:)))/(max(I(:))-min(I(:)));
            I4(Lws == 0 | imdilate(ImaxCombClean,se)) = 1;
            figure, imagesc(I4), title('Markers and object boundaries superimposed on original image (I4)')
            colormap gray;
        end
        
        % Post processing:
        % Remove those regions that don't have maxima:
        labelsWS = bwlabel(Lws>0, 4);
        if doPlot==1
            figure, imagesc(labelsWS), title('Nuclei labels')
        end
        
        maxLabelWS=max(labelsWS(:));     
        maxCombLabel=labelsWS(ImaxCombClean(:)==1);
        
        badRIOsID=find(~ismember(1:maxLabelWS,maxCombLabel));
        labelsWS(ismember(labelsWS,badRIOsID))=0;
        if doPlot==1
            figure, imagesc(labelsWS), title('Nuclei labels')
        end
        
        % Remove those regions wich are either too large (detection in the
        % noise, super large clumps) or the ones which are too small (close
        % to noise, insignificant max):
        
        % This needs more thought, instead one could apply an upper
        % threshold of direct neighbors to one regions to get the large
        % filling areas away
        applThr=1;
        if applThr==1
            display('Size thresholds have been applied!')
            labelsWS = bwlabel(labelsWS>0, 4);
            regions = regionprops(labelsWS, 'area');
            area    = vertcat((regions(:).Area));

            thrMax=(3*r)^2*pi;
            thrMin=-Inf;
            badRIOsID=vertcat(find(area>thrMax),find(area<thrMin));
            labelsWS(ismember(labelsWS,badRIOsID))=0;
        end

        % Fill holes:
        [bwBnDMaskStr, bwLabels]=bwboundaries(labelsWS>1,4,'noHoles');
        allBndPts=vertcat(bwBnDMaskStr{:});
        indBndPts=sub2ind(size(Icrop),allBndPts(:,1),allBndPts(:,2));
        bwBnDMask=zeros(size(Icrop));
        bwBnDMask(indBndPts)=1;
        
        if doPlot==1
            figure, imagesc(bwLabels), title('Nuclei labels')

            se = strel('disk', round(r/4));
            I5 = (I-min(I(:)))/(max(I(:))-min(I(:)));        
            I5(bwBnDMask == 1) = 1;
            I5(imdilate(ImaxCombClean,se)==1)=0.5;
            figure, imagesc(I5), title('Nuclei labels')
        
            I6 = (I-min(I(:)))/(max(I(:))-min(I(:)));        
            I6(bwBnDMask == 1) = 1;
            figure, imagesc(I6), title('Nuclei labels')
        end
    end  
    
    %**********************************************************************
    % 6. Step: Prepare the output                                         *
    %**********************************************************************
    
    nuclei(frame).pos=nucPos+dPix(frame);
    
    nuclei(frame).img=false(size(currImage));    
    linId=sub2ind(size(currImage),nuclei(frame).pos(:,2),nuclei(frame).pos(:,1));
    nuclei(frame).img(linId)=true;
    
    movieInfo(frame).xCoord(:,1)=nuclei(frame).pos(:,1);
    movieInfo(frame).xCoord(:,2)=r/2;
    movieInfo(frame).yCoord(:,1)=nuclei(frame).pos(:,2);
    movieInfo(frame).yCoord(:,2)=r/2;
    movieInfo(frame).amp(:,1)=currImage(linId);
    movieInfo(frame).amp(:,2)=0;
    
%     se   = strel('disk', round(r));
%     nucImgOrgDil = imdilate(nuclei(frame).img,se);
%     Idspl=currImage;
%     Idspl(nucImgOrgDil == 1) = 0;
%     figure(1), imagesc(Idspl), title('Extended nuclei');
%     colormap jet;
    
    if nargout>3
        allBndPtsOrg=allBndPts+dPix(frame);
        indBndPtsOrg=sub2ind(size(currImage),allBndPtsOrg(:,1),allBndPtsOrg(:,2));
        bwBnDMaskOrg=zeros(size(currImage));
        bwBnDMaskOrg(indBndPtsOrg)=1;
        
        [~, bwLabelsOrg]=bwboundaries(bwBnDMaskOrg>0,4,'noHoles');
        
        Idspl = (currImage-min(currImage(:)))/(max(currImage(:))-min(currImage(:)));
        Idspl(bwBnDMaskOrg==1) = 1;
        figure, imagesc(Idspl), title('Nuclei boundaries')
        
        display('Be aware of the following facts:')
        display('1.) (nucPosOrg, nucImgOrg) and (bwBnDMaskOrg, bwLabelsOrg) dont perfectly fit:');
        display('2.) detected nuclei need not to be located in boundaries!');
        display('3.) boundaries contain at least (but maybe more than) one nucleus!');
    end
end