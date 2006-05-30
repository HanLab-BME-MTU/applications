function testing_evaluate11(idlisttrack, debugData, dataProperties)
%TESTING_EVALUATE11 evaluates the results of test11
%
% SYNOPSIS: testing_evaluate11(idlisttrack, debugData, dataProperties)
%
% INPUT
%
% OUTPUT
%
% REMARKS
%
% created with MATLAB ver.: 7.2.0.232 (R2006a) on Windows_NT
%
% created by: jdorn
% DATE: 27-May-2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% prepare collecting data
positions = loadPositions(902);
positions = positions(:,[2,1,3],:);
pix2mu = [dataProperties.PIXELSIZE_XY,dataProperties.PIXELSIZE_XY,dataProperties.PIXELSIZE_Z];
nTags = size(positions,3);


id2posIdx = zeros(nTags,1);
for i=1:nTags
    id2posIdx(i) = str2double(idlisttrack(1).stats.labelcolor{i});
end


trackResults = debugData.trackResults;
numSources = size(debugData.trackResults,3);
nTimepoints = length(idlisttrack);

% sourcePos:
% t, iTag, source XYZ, subpixelposition XYZ (groundTruth), sourceDelta XYZ,
% sqrt(sourceVar XYZ)

% deltaPos
% t, iTag, tSource, delta, trueDelta, sourceSubPixel, targetSubPixel, deltaDelta,
% sqrt(deltaVar), iter, success

% trackPos
% t, iTag, nSources, pos, true pos (subpixel), deltaPos, sqrt(var)

sourcePos = zeros(nTimepoints * nTags,14);
deltaPos = zeros(nTimepoints * nTags * numSources,24);
trackPos = zeros(nTimepoints * nTags,15);
sourceCt = 0;
deltaCt = 0;
trackCt = 0;


for t=1:nTimepoints
    if ~isempty(idlisttrack(t).linklist)
        trackQDiag = diag(idlisttrack(t).info.totalQ_Pix);
        for i=1:nTags
            iTag = id2posIdx(i);

            truePos = positions(t,:,i);
            

            if isempty(debugData.trackResults(t,1,1).info) ||...
                debugData.trackResults(t,1,1).info(2) == 1;
                % it's a source
                sourceCt = sourceCt + 1;

                sourcePos(sourceCt,1) = t;
                sourcePos(sourceCt,2) = iTag;
                sourcePos(sourceCt,3:5) = debugData.trackResults(t,iTag,1).sourcePos;
                sourcePos(sourceCt,12:14) = sqrt(trackResults(t,iTag,1).sourceVar);

                sourcePos(sourceCt,9:11) = truePos - sourcePos(sourceCt,3:5);
                sourcePos(sourceCt,6:8) = truePos - floor(truePos);

                % start looking at start/end from #2
                s1=2;
            else
                s1=1;
            end

            % loop tracking
            for s = s1:numSources
                deltaCt = deltaCt + 1;
                deltaPos(deltaCt,1) = t;
                deltaPos(deltaCt,2) = iTag;

                % if clause for error correction
                if size(debugData.trackResults(t,iTag,s).startEndDelta,2) == 6
                    deltaPos(deltaCt,4:6) = debugData.trackResults(t,iTag,s).startEndDelta(2,4:6);
                else
                deltaPos(deltaCt,4:6) = debugData.trackResults(t,iTag,s).startEndDelta(2,:);
                end
                deltaPos(deltaCt,19:21) = sqrt(debugData.trackResults(t,iTag,s).deltaVar);

                sourceT = debugData.trackResults(t,1,s).info(3);
                deltaPos(deltaCt,3) = sourceT;

                truePosS = positions(sourceT,:,i);

                deltaPos(deltaCt,7:9) = truePos - truePosS;
                deltaPos(deltaCt,10:12) = truePosS - floor(truePosS);
                deltaPos(deltaCt,13:15) = truePos - floor(truePos);
                deltaPos(deltaCt,16:18) = deltaPos(deltaCt,7:9)-deltaPos(deltaCt,4:6);
                
                % remember iter,success (better at beginning, but I don't want
                % to change all the indices now
                deltaPos(deltaCt,23:24) = debugData.trackResults(t,1,s).info(6:7);
            end

            % read idlisttrack if not estimated tag
            if idlisttrack(t).linklist(iTag,3) ~= 1
                trackCt = trackCt + 1;
                trackPos(trackCt,1) = t;
                trackPos(trackCt,2) = iTag;
                trackPos(trackCt,3) = length(idlisttrack(t).info.sourceList);

                trackPos(trackCt,4:6) = idlisttrack(t).linklist(iTag,9:11)./pix2mu;
                trackPos(trackCt,10:12) = truePos - trackPos(trackCt,4:6);
                trackPos(trackCt,7:9) = truePos - floor(truePos);


                trackPos(trackCt,13:15) = trackQDiag( (iTag-1)*3+1:iTag*3 );

            end
        end
    end
end

% plot
figure

