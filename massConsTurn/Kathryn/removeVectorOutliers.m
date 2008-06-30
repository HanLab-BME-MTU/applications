function removeVectorOutliers(runInfo)
% removeVectorOutliers removes outliers from a vector field

% runInfo is a structure containing at least these fields:
%
% thresh: search radius for distance matrix between all the vector origin
% pairs in raw data set. it should be big enough so that the average
% vector has several neighbors within thresh. thresh will be optimized
% based on vector density, so that on average there will be 7 neighbors
% for a given vector.
%
% nNeighborLowerLimit: vectors with a smaller number of neighbors will be
% considered outliers and removed from the filtered set. 3 is recommended.
%
% dirDev: mean directional deviation allowed (in degrees) for a vector and
% its neighbors. we find the mean(cos(theta)) for the vector and its
% neighbors and compare it to the corresponding value calculated from
% dirDev. vectors which deviate more than this value are considered
% outliers and removed from the filtered set.


% check input and assign parameters
if nargin<1
    error('removeVectorOutliers: Not enough input parameters')
end
if ~isstruct(runInfo)
    runInfo=struct;
end

if ~isfield(runInfo,'imDir') || ~isfield(runInfo,'anDir')
    error('removeVectorOutliers: runInfo should contain fields imDir and anDir');
else
    [runInfo.anDir] = formatPath(runInfo.anDir);
    [runInfo.imDir] = formatPath(runInfo.imDir);
end

if isfield(runInfo,'thresh')
    thresh=runInfo.thresh;
else
    thresh=30;
end
if isfield(runInfo,'nNeighborLowerLimit')
    nNeighborLowerLimit=runInfo.nNeighborLowerLimit;
else
    nNeighborLowerLimit=3;
end
if isfield(runInfo,'dirDev')
    dirDev=runInfo.dirDev;
else
    dirDev=180;
end

% flow tracking directory
corrDir=[runInfo.anDir filesep 'corr'];
if ~isdir(corrDir)
    mkdir(corrDir);
end

vecDir=[corrDir filesep 'filt'];
if isdir(vecDir)
    rmdir(vecDir,'s');
end
mkdir(vecDir);
mkdir([vecDir filesep 'vecOutliers']);
mkdir([vecDir filesep 'vectorCoverageMask']);

% image negative directory
negImDir=[runInfo.imDir filesep 'norm' filesep 'negEdgeTifs'];
if ~isdir(negImDir)
    mkdir(negImDir);
end
listOfImages=searchFiles('neg_edge',[],negImDir,0);
nImages=size(listOfImages,1);


% check for flow track format
% searchFiles not case sensitive; if you change this, be careful that names conform to
% standards

% vecPxyVxyAllFms: nFlwTrcks x 1 cell containing the tail positions
% and x- and y-components of all vectors
listFilesLin = searchFiles('flowTrack(\w*)_(\w*).mat',[],corrDir,0);
listFilesCyrus = searchFiles('(\w*)FlowTrack.mat','_',corrDir,0);
listFilesJames= searchFiles('vectors(\w*).mat',[],corrDir,0);

if ~isempty(listFilesLin(:,1)) && isempty(listFilesCyrus(:,1)) % Lin format
    runInfo.nFlwTrcks=length(listFilesLin(:,1)); %number of flowtracks

    % fill cell with measured velocity info from all tracks
    vecPxyVxyAllFms=cell(runInfo.nFlwTrcks,1);
    for n=1:runInfo.nFlwTrcks
        fileName=[char(listFilesLin(n,2)) filesep char(listFilesLin(n,1))];
        flwTrck=load(fileName);
        vecPxyVxyAllFms{n}=[flwTrck.flowTrack.p{1,1} flwTrck.flowTrack.v{1,1}];
    end

elseif isempty(listFilesLin(:,1)) && ~isempty(listFilesCyrus(:,1)) % Cyrus format
    fileName=[char(listFilesCyrus(1,2)) filesep char(listFilesCyrus(1,1))];
    flwTrck=load(fileName);
    flwTrck=flwTrck.flowHistory;
    runInfo.nFlwTrcks=length(flwTrck); %number of flowtracks

    % fill cell with measured velocity info from all tracks
    vecPxyVxyAllFms=cell(runInfo.nFlwTrcks,1);
    for n=1:runInfo.nFlwTrcks
        vecPxyVxyAllFms{n}=[flwTrck(1,n).p flwTrck(1,n).v];
    end

elseif ~isempty(listFilesJames(:,1))
    runInfo.nFlwTrcks=length(listFilesJames(:,1)); %number of flowtracks

    % fill cell with measured velocity info from all tracks
    vecPxyVxyAllFms=cell(runInfo.nFlwTrcks,1);

    for n=1:runInfo.nFlwTrcks
        fileName=[char(listFilesJames(n,2)) filesep char(listFilesJames(n,1))];
        flwTrck=load(fileName);
        vecPxyVxyAllFms{n}=[flwTrck.vectors(:,2:-1:1) flwTrck.vectors(:,4:-1:3)-flwTrck.vectors(:,2:-1:1)];
    end

else
    disp('POLYDEPOLYMAP: flow tracks missing from /corr or have unknown format')
end



nFlwTrcks=runInfo.nFlwTrcks;
filteredFlow(nFlwTrcks)=struct('pyxRaw',[],'vyxRaw',[],'pyxFilt',[],'vyxFilt',[]);
s=length(num2str(nFlwTrcks));
strg=sprintf('%%.%dd',s);

runInfo.nImages=min(nImages,nFlwTrcks);

h=get(0,'CurrentFigure');
if h==1
    close(h)
end
figure(1);
for i=1:runInfo.nImages

    % extract the vectors for frm1, separate into known vector tail
    % positions and y and x components
    vecPxyVxy=vecPxyVxyAllFms{i,1};
    pyxVecKnown=vecPxyVxy(:,2:-1:1);
    vyxVecKnown=vecPxyVxy(:,4:-1:3);
    pyxVecKnown(isnan(vyxVecKnown(:,1)),:)=[];
    vyxVecKnown(isnan(vyxVecKnown(:,1)),:)=[];


    % get the distance between every vector pair in pyxVecKnown
    D=createDistanceMatrix(pyxVecKnown,pyxVecKnown);

    % inRange is 1 where vector (in y-direction) fits criteria
    inRange=(D>0 & D<=thresh);
    nNeighbors=sum(inRange,2);

    % adapt thresh value until each vector has an average of 7 neighbors within
    % thresh
    meanNumNeighbors=mean(nNeighbors);
    while meanNumNeighbors>7
        thresh=thresh-1;
        inRange=D>0 & D<=thresh;
        nNeighbors=sum(inRange,2);
        meanNumNeighbors=mean(nNeighbors);
    end

    % remove vectors without enough neighbors
    inRange(nNeighbors<nNeighborLowerLimit,:)=0;
    inRange=swapMaskValues(inRange,0,nan); % use NaN's instead of 0's

    % using dot(A,B)=|A||B|cos(theta), we can find how coherent the
    % vectors are.  coherency is measured by cos(theta) for each pair
    mag=sqrt(sum(vyxVecKnown.^2,2)); % vector magnitudes
    AB=mag*mag'; % matrix containing |A||B|
    My=vyxVecKnown(:,1)*vyxVecKnown(:,1)'; % product of y-components
    Mx=vyxVecKnown(:,2)*vyxVecKnown(:,2)'; % product of x-components
    cosTheta=(My+Mx)./AB;

    % convert dirDev to radians and take cosine
    cutOff=cos(dirDev*(pi()/180));

    % remove incoherent vectors
    inRange(nanmean(cosTheta.*inRange,2)<cutOff,:)=nan;

    % retain any vector that still has a neighbor in inRange
    idxFilt=nansum(inRange,2)>=1;
    pyxFiltered=pyxVecKnown(idxFilt,:);
    vyxFiltered=vyxVecKnown(idxFilt,:);

    runInfo.meanNN=round(nanmean(nanmin(D.*inRange,[],2)));

    % load the image negative
    imgName=[char(listOfImages(i,2)) filesep char(listOfImages(i,1))];
    img=double(imread(imgName));
    [runInfo.imL, runInfo.imW]=size(img);

    % make mask with 1 on all pixels where a vector tail is located, then
    % dilate/erode
    vectorCoverageMask=zeros(runInfo.imL,runInfo.imW);
    vectorCoverageMask(xy2index(pyxFiltered(:,2),pyxFiltered(:,1),runInfo.imL,runInfo.imW))=1;
    maxNN=ceil(nanmax(nanmin(D.*inRange,[],2)));

    vectorCoverageMask = imdilate(vectorCoverageMask,strel('disk',3*maxNN));
    vectorCoverageMask = imerode(vectorCoverageMask,strel('disk',maxNN));
    vectorCoverageMask=logical(vectorCoverageMask);


    edgePix=find(bwmorph(vectorCoverageMask,'remove'));
    img(edgePix)=round(max(img(:))/2); % gray border around coverage mask


    imshow(img,[]);
    hold on
    quiver(pyxVecKnown(:,2),pyxVecKnown(:,1),vyxVecKnown(:,2),vyxVecKnown(:,1),0,'r');
    quiver(pyxFiltered(:,2),pyxFiltered(:,1),vyxFiltered(:,2),vyxFiltered(:,1),0,'b');
    drawnow
    indxStr1=sprintf(strg,i);
    h=get(0,'CurrentFigure');
    saveas(h,[vecDir filesep 'vecOutliers' filesep 'overlay' indxStr1],'tif');
    close

    filteredFlow(i).pyxRaw=pyxVecKnown;
    filteredFlow(i).vyxRaw=vyxVecKnown;
    filteredFlow(i).pyxFilt=pyxFiltered;
    filteredFlow(i).vyxFilt=vyxFiltered;

    imwrite(vectorCoverageMask,[vecDir filesep 'vectorCoverageMask' filesep 'vectorCoverageMask' indxStr1 '.tif']);
    

end


save([vecDir filesep 'filteredFlow'],'filteredFlow');
save([runInfo.anDir filesep 'runInfo'],'runInfo');
