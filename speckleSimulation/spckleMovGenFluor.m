function [params,py1,px1] = spckleMovGenFluor(params,not2keep)
% generate fluorophores in area of image and border for t=0

switch params.nModel
    case {1, 2, 3}
    % for actin, the minimum repeat of the double-stranded helix is 35.75 nm and contains
    % 13 monomers (Biophys J. 2003 July; 85(1): 27–35)
    areaDensNM=params.areaDens*0.001*(13/35.75); % f-actin monomer density (subunits/nm^2)
    avgActinPerSquare=areaDensNM*params.sampleScale^2; % avg number subunits per sampling box
    avgLabeledActinPerSquare=params.percentLabeledActin*avgActinPerSquare; % avg number labeled subunits per box

    % calculate how many fluorophores to use initially
    nFluor=ceil(avgLabeledActinPerSquare*not2keep.nSampleSquares);

    py1=not2keep.nSamplesL*params.sampleScale*rand(nFluor,1); % y-coordinate of each fluorophore
    px1=not2keep.nSamplesW*params.sampleScale*rand(nFluor,1); % x-coordinate of each fluorophore

    case 4
    % convert speeds to nm/s
    fastMean = params.fastFlowMean*1000/60;
    fastSigma = params.fastFlowSigma*1000/60;
    slowMean = params.slowSpeedMean*1000/60;
    slowSigma = params.slowSpeedSigma*1000/60;

    % num MTs and spaces in each bundle/space/singles/space repeat
    nInRepeat = params.nPerBundle + params.nSpacers1 + params.nSingles + params.nSpacers2;

    nPossiblePos = floor((params.imL*params.pixNM)/24); % num of parallel MTs that can fit in the image (horizontal)
    yPositions = 24*[1:nPossiblePos]'-12; % y-coordinates of fluorophore placed in center of the 24nm-wide MT

    % tag =
    %  -1        empty spaces
    %   0        single MTs
    %   1,2,3... bundles
    %   each MT in a bundle gets the same tag

    tag = floor([0:nPossiblePos-1]'./nInRepeat)+1;

    for i=1:params.nSpacers1 % empty spaces get tag of -1
        tag(params.nPerBundle+i:nInRepeat:end)=-1;
    end

    for i=1:params.nSingles % single MTs get tag of 0
        tag(params.nPerBundle+params.nSpacers1+i:nInRepeat:end)=0;
    end

    for i=1:params.nSpacers2 % empty spaces get tag of -1
        tag(params.nPerBundle+params.nSpacers1+params.nSingles+i:nInRepeat:end)=-1;
    end

    fluorData.tag=tag; % keep the vector of tags

    nBundles=sum(unique(tag)>0); % count number of bundles
    bundleTags=unique(tag(tag>0)); % get a list of the tags used to identify the bundles

    nIndivMTs=sum(tag==0);  % count number of individual MTs
    indivTags=find(tag==0); % get tags used to identify the individual MTs

    labelingRatio=zeros(nPossiblePos,1); % initialize labeling ratio vector
    labelingRatio(tag~=0 & tag~=-1)=params.effectiveLabeling; % bundle MTs get effective labeling value
    labelingRatio(indivTags)=params.fractionTuLabeled; % individual MTs get lower value

    flowDir=zeros(nPossiblePos,1); % initialize flow direction vector
    flowVel=zeros(nPossiblePos,1); % initialize flow speed vector - will choose from distribution

    % get speed and direction for MTs in each bundle
    % in a given bundle, all fluorophores move at constant speed (chosen from slow
    % distribution) in the same direction (chosen randomly)
    for i=1:nBundles
        flowVel(tag==bundleTags(i))=slowMean+slowSigma*randn;
        if rand<=0.5 % equal chance for motion in either direction
            flowDir(tag==bundleTags(i))=-1;
        else
            flowDir(tag==bundleTags(i))=1;
        end
    end
    % get speed and direction for individual MTs (most in slow population)
    nFast = round(params.fractionFastSingleMTs*nIndivMTs);
    randList=randperm(nIndivMTs);
    fastMTIndices=indivTags(randList(1:nFast));
    slowMTIndices=setdiff(indivTags,fastMTIndices);

    for i=1:nFast
        flowVel(fastMTIndices(i))=fastMean+fastSigma*randn;
        if rand<=0.5 %assign direction
            flowDir(fastMTIndices(i))=-1;
        else
            flowDir(fastMTIndices(i))=1;
        end
    end
    for i=1:length(slowMTIndices)
        flowVel(slowMTIndices(i))=slowMean+slowSigma*randn;
        if rand<=0.5 %assign direction
            flowDir(slowMTIndices(i))=-1;
        else
            flowDir(slowMTIndices(i))=1;
        end
    end

    % velocity is product of speed and direction vector
    flowVel=flowVel.*flowDir;

    nmPerSecFlowSpeed=max(abs(flowVel)); % how many nm the fastest MT moves per second
    not2keep.nmPerFrame=nmPerSecFlowSpeed*params.nSecPerFrame; % how many nm the MT flows per frame of the movie

    nTuSubunitsPerRow=round((not2keep.imgWnm/8)); % number of tubulin subunits per MT in image
    nToLabel=round(nTuSubunitsPerRow*labelingRatio); % vector of number of labeled tu subunits to have per MT

    px1=not2keep.imgWnm*rand(sum(nToLabel),1); % positions of all fluorophores along width in nm
    py1=zeros(length(px1),1); % initialize y positions
    params.mod4vx=zeros(length(px1),1); % initialize velocity vector
    
    for i=1:length(nToLabel)
        if i==1
            counter=1;
        else
            counter=sum(nToLabel(1:i-1))+1;
        end
        py1(counter:counter-1+nToLabel(i))=yPositions(i);
        params.mod4vx(counter:counter-1+nToLabel(i))=flowVel(i);
    end
    
    otherwise
    % future models
end

