function anaDat = adgui_calc_anaDat(idlist,dataProperties,lastResult)
%ADGUI_CALC_ANADAT brings the data stored in the idlist into a form that is convenient for subsequent data analysis.
%
%SYNOPSIS anaDat = adgui_calc_anaDat(idlist,dataProperties,lastResult)
%
%INPUT idlist/dataProperties  data structures generated with the analyzeMoviesGUI/runCtBatch-package
%      lastResult             indicates whether the data structure has run through the tracking module yet 
%                               (in that case, lastResult is 'idlisttrack')
%
%OUTPUT anaDat: Structure containing the fields
%         (1).info
%               .name:       name of project
%               .idlisttype  =lastResult (either 'idlist' or 'idlisttrack')
%               .nTags:      number of tags in movie
%               .labelColor  names of tags 1:n (cell array)
%               .created     date created (nowString)
%               .dataProperties dataProperties of project
%
%       (1:n).coord                  tag coordinates
%       (1:n).displacementVectorsN   nTags by 3 list of normed displacement vectors
%       (1:n).displacement           norms of the above vectors
%       (1:n).centroid               centroid of the tags
%       (1:n).distanceVectorMatrixN  nTags by nTags by 3 distance vector matrix
%       (1:n).distanceMatrix         nTags by nTags distanceMatrix
%       (1:n).time                   mean time of frame
%       (1:n).sigmaTime              (max-min)/2 time of frame
%       (1:n).timePoint              timePoint number (entry in idlist)
%
%       (1:n).stats
%               .qMatrix        nTags+1 by nTags+1 blocks of 3 by 3
%                                   Q-matrices (1/sigma(ik)^2)
%               .noise          noise estimation for each tag
%
%based on analyzeIdlist by dT
%
%c: Jonas, 05/03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%init counter
gT = 1;

%calculate time
acqTime = mean(dataProperties.frameTime,2);
sigmaTime = (max(dataProperties.frameTime,[],2)-min(dataProperties.frameTime,[],2))/2;
%get nTimePoints
nTimePoints = size(idlist,2);

%find first nonempty entry in idlist
startT=1;
while(isempty(idlist(startT).linklist) & startT<nTimePoints)
    startT=startT+1;
end;
if startT > nTimePoints
    error('empty idlist')
end


%write idlist information
anaDat(1).info.name = dataProperties.name;
anaDat(1).info.idlisttype = lastResult;
anaDat(1).info.nTags = size(idlist(startT).linklist,1);
anaDat(1).info.labelColor = idlist(1).stats.labelcolor(1:anaDat(1).info.nTags); %take no risks with old idlists
anaDat(1).info.created = nowString;
anaDat(1).info.dataProperties = dataProperties;

%calculate transformation matrix pix->mu (for q-matrices)
pix2muMat = [];
p2mM = diag([dataProperties.PIXELSIZE_XY,dataProperties.PIXELSIZE_XY,dataProperties.PIXELSIZE_Z]);
for i=1:anaDat(1).info.nTags
    pix2muMat = blkdiag(pix2muMat,p2mM);
end

s  =  warning('query', 'all');
warning off MATLAB:divideByZero

%fill in first entry in anaDat (gT=1)
idlist(startT).linklist = sortrows(idlist(startT).linklist,4); %sort idlist 
anaDat(gT).coord = idlist(startT).linklist(:,9:11); %tag coordinates
anaDat(gT).displacementVectorsN = []; %displacement of tags 1->2 will be stored in 1 later
anaDat(gT).displacement = []; %same here
anaDat(gT).centroid = idlist(startT).centroid; %centroid
[anaDat(gT).distanceMatrix,distanceVectorMatrix] =...
    distMat(anaDat(gT).coord); %nTags by nTags matrix of distance between points 
anaDat(gT).distanceVectorMatrixN = ... %nTags by nTags by 3 matrix of normed distance vectors
    distanceVectorMatrix./cat(3,anaDat(gT).distanceMatrix,anaDat(gT).distanceMatrix,anaDat(gT).distanceMatrix);
anaDat(gT).time = acqTime(startT); %true acquisition time
anaDat(gT).timePoint = startT; 

%info and stats are the other way round in idlist, but this is more logical.
%QMatrices are already calculated in microns, as for all subsequent analysis we
%calculate in that space (and transformation of distance is non-linear!)
if isfield(idlist(startT).info,'trackQ_Pix') & ~isempty(idlist(startT).info.trackQ_Pix)
    %tracker gives only relative tracking error. do squared addition of
    %source error with relative tracking error -> will crash if multiple
    %sources are being used
    sourceT = str2double(idlist(startT).info.trackerMessage.source);
    anaDat(gT).stats.qMatrix = pix2muMat*(idlist(startT).info.trackQ_Pix+idlist(sourceT).info.detectQ_Pix)*pix2muMat; 
else
    anaDat(gT).stats.qMatrix = pix2muMat*idlist(startT).info.detectQ_Pix*pix2muMat;
end
%add QMatrix of CoM
allQ = zeros(3,3,anaDat(1).info.nTags);
for i = 1:anaDat(1).info.nTags %read out q-submatrices
    allQ(:,:,i) = anaDat(gT).stats.qMatrix( (i-1)*3+1:i*3, (i-1)*3+1:i*3 );
end
%take the mean of all of them
comQ = mean(allQ,3)/anaDat(1).info.nTags^2;
anaDat(gT).stats.qMatrix = blkdiag(anaDat(gT).stats.qMatrix,comQ);

%the old noise is likely to be not correct - but we want to keep the
%program working - and there's the simulation
if size(idlist(startT).linklist,2)>11
    anaDat(gT).stats.noise = idlist(startT).linklist(:,12)';
    noiseInLL = 1;
else
    anaDat(gT).stats.noise = idlist(startT).info.noise;
    noiseInLL = 0;
end
anaDat(gT).sigmaTime = sigmaTime(startT);


for t=startT+1:nTimePoints
    if ~isempty(idlist(t).linklist)
        idlist(t).linklist = sortrows(idlist(t).linklist,4); %sort idlist 
        
        %get displacement (t1->t2 is stored in t1)
        displacementVectors = idlist(t).linklist(:,9:11)-anaDat(gT).coord;
        [anaDat(gT).displacement, anaDat(gT).displacementVectorsN] = normList(displacementVectors);
        
        %fill in values for t2
        gT=gT+1;
        
        anaDat(gT).coord = idlist(t).linklist(:,9:11);
        anaDat(gT).centroid = idlist(t).centroid;
        [anaDat(gT).distanceMatrix,distanceVectorMatrix] =...
            distMat(anaDat(gT).coord); %nTags by nTags matrix of distance between points 
        anaDat(gT).distanceVectorMatrixN = ... %nTags by nTags by 3 matrix of normed distance vectors
            distanceVectorMatrix./cat(3,anaDat(gT).distanceMatrix,anaDat(gT).distanceMatrix,anaDat(gT).distanceMatrix);
        anaDat(gT).time = acqTime(t);
        anaDat(gT).timePoint = t;
        %get QMatrix and transform to microns
        if isfield(idlist(t).info,'trackQ_Pix') & ~isempty(idlist(t).info.trackQ_Pix)
            sourceT = str2double(idlist(t).info.trackerMessage.source);
            anaDat(gT).stats.qMatrix = pix2muMat*(idlist(t).info.trackQ_Pix+idlist(sourceT).info.detectQ_Pix)*pix2muMat;  
        else
            anaDat(gT).stats.qMatrix = pix2muMat*idlist(t).info.detectQ_Pix*pix2muMat;
        end
        %QMatrix for centroid: mean of all q-matrices divided by n^2
        allQ = zeros(3,3,anaDat(1).info.nTags);
        for i = 1:anaDat(1).info.nTags %read out q-submatrices
            allQ(:,:,i) = anaDat(gT).stats.qMatrix( (i-1)*3+1:i*3, (i-1)*3+1:i*3 );
        end
        comQ = mean(allQ,3)/anaDat(1).info.nTags^2;
        anaDat(gT).stats.qMatrix = blkdiag(anaDat(gT).stats.qMatrix,comQ);
        
        if noiseInLL
            anaDat(gT).stats.noise = idlist(t).linklist(:,12)';
        else
            anaDat(gT).stats.noise = idlist(t).info.noise;
        end
        anaDat(gT).sigmaTime = sigmaTime(t);
    end;
end;

warning(s); %turn back on the warnings

