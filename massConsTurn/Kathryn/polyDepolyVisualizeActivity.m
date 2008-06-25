function polyDepolyVisualizeActivity(runInfo,mapDir,runInfoMinMaxValues,cMap,gamma,mkMov)
% POLYDEPOLYVISUALIZEACTIVITY makes a vector overlay movie of polyDepoly
%
% SYNOPSIS: polyDepolyVisualizeActivity(runInfo,mapDir,runInfoMinMaxValues,cMap,gamma,mkMov)
%
% INPUT:
%    runInfo             : structure created by turnoverMap
%    mapDir              : where poly depoly maps are stored
%    runInfoMinMaxValues : 1 if you want to use the min/max values for
%                          poly-depoly recorded in runInfo; 0 if you want
%                          to calculate the values from scratch (e.g. time
%                          averaged ones will have a smaller range than the
%                          individual frame stack, so using 0 maximizes
%                          their colors)
%    cMap                : 'iso' if isomorphic red-green, 'jet' for rainbow
%    gamma               : value for gamma correction (<1 boosts low signals)
%                          (default = 1: no change)
%
% OUTPUT:
%    a new directory is created (\turn\interp\mapTifs) where output tifs
%    and movie are stored. several parameters are saved in the .mat file
%    visualizationParams
%
% underlying normalized images are brightened using imAdjust for purposes
% of movie making.
%
% MATLAB VERSION (originally written on): 7.2.0.232 (R2006a) Windows XP
%
% USERNAME: kathomps
% DATE: 12-Apr-2006
%

warningState=warning;
warning('off','MATLAB:intConvertNonIntVal')

% retrieve raw vectors in cyan, filtered set in yellow (only outliers show
% up as cyan)
filteredFlow=load([runInfo.turnDir filesep 'filteredFlow.mat']);
filteredFlow=filteredFlow.filteredFlow;

if isequal(cMap,'iso')
    % create the color map with black at the end
    % this function makes a perceptual color map
    clrMap=isomorphicColormap('g/r',64);
    clrMap=[clrMap; [0 0 0]];
    opacity=.5;
elseif isequal(cMap,'jet')
    clrMap=colormap('jet');
    opacity=.5;
else
    error('POLYDEPOLYVIZUALIZEACTIVITY: colormap not supported')
end

% get list of activity maps and count them
[listOfMaps] = searchFiles('polyDepoly','tif',mapDir,0);
nMaps=size(listOfMaps,1);
s=length(num2str(nMaps));
strg=sprintf('%%.%dd',s);

% frameIdx gives which frames to use for the underlying images
% the function first looks to see if the maps are in a time-averaged
% directory and loads the parameters used to generate it
if exist([mapDir filesep 'tmAvgParams.mat'])
    load([mapDir filesep 'tmAvgParams.mat'])
    firstFrm=[startFrm:timeStepSize:endFrm]; lastFrm=firstFrm+nFrms2Avg-1;
    firstFrm(lastFrm>endFrm)=[];             lastFrm(lastFrm>endFrm)=[];
    frameIdx=firstFrm+floor(nFrms2Avg/2);
else % not time-averaged directory; use sequential images
    frameIdx=1:nMaps;
end

% get list of normalized images
[listOfImages] = searchFiles('norm',[],[runInfo.anDir filesep 'norm' filesep 'normMats'],0);

% get list of cell masks
cmDir=[runInfo.anDir filesep 'edge' filesep 'cell_mask'];
[listOfCellMasks] = searchFiles('.tif',[],cmDir,0);

% get min/max over whole image series if not using global values from batch
if runInfoMinMaxValues==0
    minValue=0;
    maxValue=0;
    for i=1:nMaps
        % load the image
        iMap=load([char(listOfMaps(i,2)) filesep char(listOfMaps(i,1))]);
        polyDepoly=eval(['iMap.' char(fieldnames(iMap))]);

        minValue=min(minValue,nanmin(polyDepoly(:)));
        maxValue=max(maxValue,nanmax(polyDepoly(:)));
    end
else
    minValue=runInfo.polyDepolyMovieMin;
    maxValue=runInfo.polyDepolyMovieMax;
end

% m is the upper bound for normalization to color map range
m=max(abs([minValue maxValue]));

% if gamma correction, apply it to the min/max values so scaling is right
if gamma~=1
    m=m.^gamma;
end

% make output directory
outDir=[mapDir filesep 'mapTifs_' num2str(m) '_' num2str(gamma)];
if ~isdir(outDir)
    mkdir(outDir);
else
    delete([outDir filesep '*tif'])
end


for i=1:nMaps
    % load the maps - these can be read sequentially b/c they're already in order (don't need frameIdx)
    iMap=load([char(listOfMaps(i,2)) filesep char(listOfMaps(i,1))]);
    polyDepoly=eval(['iMap.' char(fieldnames(iMap))]);
    
    
    % gamma-correct pixel values here (boost low values if gamma < 1)
    if gamma~=1
        sign=zeros(size(polyDepoly));
        sign(polyDepoly<0)=-1;
        sign(polyDepoly>=0)=1;
        polyDepoly=sign.*(abs(polyDepoly).^gamma);
    end
    
    % load cell mask (need frameIdx)
    cMask=double(imread([char(listOfCellMasks(frameIdx(i),2)) filesep char(listOfCellMasks(i,1))]));
    cMask=swapMaskValues(cMask,0,nan);
    polyDepolyMasked=cMask.*polyDepoly;

    % load the normalized image of interest (need frameIdx)
    iImg=load([char(listOfImages(frameIdx(i),2)) filesep char(listOfImages(frameIdx(i),1))]);
    im=eval(['iImg.' char(fieldnames(iImg))]);
    im=imadjust(im,[0 1],[0 1],.33);
    
    % overlay
    [img3C,pixClasses]=imDataMapOverlay(im,polyDepolyMasked,[-m,m],clrMap,opacity);

    % show edge from corresponding image in yellow (middle if time-averaged)
    img3C(runInfo.edgePix{frameIdx(i)})=1;
    img3C(runInfo.edgePix{frameIdx(i)}+runInfo.imL*runInfo.imW)=1;

    hmain=get(0,'CurrentFigure');
    if i==1 && ~isempty(hmain)
        figure(hmain+1);
    else
        figure(1);
    end
    
    % show overlay and vectors from corresponding frame (middle if time-averaged)
    imshow(img3C);
    hold on
    quiver(filteredFlow(1,frameIdx(i)).pyxRaw(:,2),filteredFlow(1,frameIdx(i)).pyxRaw(:,1),filteredFlow(1,frameIdx(i)).vyxRaw(:,2),filteredFlow(1,frameIdx(i)).vyxRaw(:,1),0,'c')
    quiver(filteredFlow(1,frameIdx(i)).pyxFilt(:,2),filteredFlow(1,frameIdx(i)).pyxFilt(:,1),filteredFlow(1,frameIdx(i)).vyxFilt(:,2),filteredFlow(1,frameIdx(i)).vyxFilt(:,1),0,'y')

    % save frames into movie
    if mkMov==1
        if i==1
            eval(['MakeQTMovie start ' [outDir filesep 'polyDepoly.mov']])
        end

        MakeQTMovie addaxes
        MakeQTMovie('framerate', 2);
        
        if i==nMaps
            MakeQTMovie finish
        end
    end

    % write the image overlay (without vectors)
    indxStr=sprintf(strg,i);
    imwrite(img3C,[outDir filesep 'polyDepoly' indxStr '.tif']);

end
warning(warningState);

save([outDir filesep 'visualizationParams'],'minValue','maxValue','m','gamma','cMap','runInfoMinMaxValues');

% for i=1:nMaps
%     % load mask and bound by fieldGeom-maxProjectionEdge
%     if isempty(isTif{1,1}) % if masks are stored as mat
%         iMask=load([char(listOfCellMasks(i,2)) filesep char(listOfCellMasks(i,1))]);
%         iMask=eval(['iMask.' char(fieldnames(iMask))]);
%     else % must be a tif
%         iMask=double(imread([char(listOfCellMasks(i,2)) filesep char(listOfCellMasks(i,1))]));
%     end
%     % bound by roi and **make mask binary**
%     boundMask=logical(iMask).*roiMask;
%     % get rid of zeros so these don't become the minValue
%     boundMask=swapMaskValues(boundMask,0,nan);
%     % load the image
%     iMap=load([char(listOfMaps(i,2)) filesep char(listOfMaps(i,1))]);
%     im=eval(['iMap.' char(fieldnames(iMap))]);
%     % backfill image with NaNs
%     im=im.*boundMask;
%
%     %make cell mask RGB size
%     cellMask=repmat(iMask,[1 1 3]);
%
%     % gamma-correct pixel values here (boost low values if gamma < 1)
%     sign=zeros(size(im));
%     sign(im<0)=-1;
%     sign(im>=0)=1;
%     im=sign.*(abs(im).^gamma);
%
%     % rescale the image
%     im(im<-m)=-m;
%     im(im>m)=m;
%     im=im+m;                    % range is 0-2m
%     im=im./(2*m);               % range is 0-1
%     im=im.*(cMapLength-1);      % range is 0-63 (if cMapLength is 64)
%     im=round(im);               % round to integers
%     im=im+1;                    % range is 1-64 (-32 to 32)
%     im(isnan(im))=cMapLength+1;
%
%     % initialize RGB arrays
%     R=zeros(runInfo.imL,runInfo.imW);
%     G=zeros(runInfo.imL,runInfo.imW);
%     B=zeros(runInfo.imL,runInfo.imW);
%
%     % fill RGB arrays with image values
%     R(:)=colorMap(im(:),1);
%     G(:)=colorMap(im(:),2);
%     B(:)=colorMap(im(:),3);
%
%     % make the cell edge white
%     R(runInfo.edgePix{i})=1;
%     G(runInfo.edgePix{i})=1;
%     B(runInfo.edgePix{i})=1;
%
%     % make the final RGB image
%     RGmap=zeros(runInfo.imL,runInfo.imW,3);
%     RGmap(:,:,1)=R;
%     RGmap(:,:,2)=G;
%     RGmap(:,:,3)=B;
%
%     % apply the cell mask
%     %RGmap=RGmap.*cellMask;
%
%     % write the image
%     indxStr=sprintf(strg,i);
%     imwrite(RGmap,[outDir filesep 'polyDepoly' indxStr '.tif']);
%     % disp(['Saving red-green tiffs: frame ' num2str(i)])
%
% end


