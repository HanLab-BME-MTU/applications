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

if nargin<1
    error('polyDepolyVisualizeActivity: Not enough input parameters')
end
if ~isstruct(runInfo)
    runInfo=struct;
end

if ~isfield(runInfo,'imDir') || ~isfield(runInfo,'anDir')
    error('polyDepolyVisualizeActivity: runInfo should contain fields imDir and anDir');
else
    [runInfo.anDir] = formatPath(runInfo.anDir);
    [runInfo.imDir] = formatPath(runInfo.imDir);
end


% flow tracking directory, retrieve flow data
filtDir=[runInfo.anDir filesep 'corr' filesep 'filt'];
temp=load([filtDir filesep 'filteredFlow.mat']);
filteredFlow=temp.filteredFlow;



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
[listOfImages] = searchFiles('norm',[],[runInfo.imDir filesep 'norm' filesep 'normMats'],0);

% get list of cell masks
cmDir=[runInfo.anDir filesep 'edge' filesep 'cell_mask'];
[listOfCellMasks] = searchFiles('.tif',[],cmDir,0);

% pixels on cell boundary
temp=load([runInfo.anDir filesep 'edgePix.mat']);
edgePix=temp.edgePix;


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

h=get(0,'CurrentFigure');
if h==1
    close(h)
end
figure(1);

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
    img3C(edgePix{frameIdx(i)})=1;
    img3C(edgePix{frameIdx(i)}+runInfo.imL*runInfo.imW)=1;


    % show overlay and vectors from corresponding frame (middle if
    % time-averaged)

    imshow(img3C);
    hold on
    quiver(filteredFlow(1,frameIdx(i)).pyxRaw(:,2),filteredFlow(1,frameIdx(i)).pyxRaw(:,1),filteredFlow(1,frameIdx(i)).vyxRaw(:,2),filteredFlow(1,frameIdx(i)).vyxRaw(:,1),0,'c')
    quiver(filteredFlow(1,frameIdx(i)).pyxFilt(:,2),filteredFlow(1,frameIdx(i)).pyxFilt(:,1),filteredFlow(1,frameIdx(i)).vyxFilt(:,2),filteredFlow(1,frameIdx(i)).vyxFilt(:,1),0,'y')
    drawnow
    indxStr=sprintf(strg,i);
    h=get(0,'CurrentFigure');
    saveas(h,[outDir filesep 'polyDepoly' indxStr],'tif');


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

    

end
close

warning(warningState);

save([outDir filesep 'visualizationParams'],'minValue','maxValue','m','gamma','cMap','runInfoMinMaxValues');