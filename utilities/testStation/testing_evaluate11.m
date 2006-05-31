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
positions = positions(:,[1,2,3],:);
pix2mu = [dataProperties.PIXELSIZE_XY,dataProperties.PIXELSIZE_XY,dataProperties.PIXELSIZE_Z];
nTags = size(positions,3);


id2posIdx = zeros(nTags,1);
for i=1:nTags
    id2posIdx(i) = str2double(idlisttrack(1).stats.labelcolor{i});
end


numSources = size(debugData.trackResults,3);
nTimepoints = length(idlisttrack);

% sourcePos:
% t, iTag, source XYZ, subpixelposition XYZ (groundTruth), sourceDelta XYZ,
% sqrt(sourceVar XYZ)

% deltaPos
% t, iTag, tSource, delta, trueDelta, sourceSubPixel, targetSubPixel, deltaDelta,
% sqrt(deltaVar), iter, success, sigma0_end

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
        for iTag=1:nTags
            iPos = id2posIdx(iTag);

            truePos = positions(t,:,iPos);
            

            if isempty(debugData.trackResults(t,1,1).info) ||...
                debugData.trackResults(t,1,1).info(2) == 1;
                % it's a source
                sourceCt = sourceCt + 1;

                sourcePos(sourceCt,1) = t;
                sourcePos(sourceCt,2) = iPos;
                sourcePos(sourceCt,3:5) = debugData.trackResults(t,iTag,1).sourcePos;
                sourcePos(sourceCt,12:14) = sqrt(debugData.trackResults(t,iTag,1).sourceVar);

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
                deltaPos(deltaCt,2) = iPos;

                % if clause for error correction
                if size(debugData.trackResults(t,iTag,s).startEndDelta,2) == 6
                    deltaPos(deltaCt,4:6) = debugData.trackResults(t,iTag,s).startEndDelta(2,4:6);
                else
                deltaPos(deltaCt,4:6) = debugData.trackResults(t,iTag,s).startEndDelta(2,:);
                end
                deltaPos(deltaCt,19:21) = sqrt(debugData.trackResults(t,iTag,s).deltaVar);

                sourceT = debugData.trackResults(t,1,s).info(3);
                deltaPos(deltaCt,3) = sourceT;

                truePosS = positions(sourceT,:,iPos);

                deltaPos(deltaCt,7:9) = truePos - truePosS;
                deltaPos(deltaCt,10:12) = truePosS - floor(truePosS);
                deltaPos(deltaCt,13:15) = truePos - floor(truePos);
                deltaPos(deltaCt,16:18) = deltaPos(deltaCt,7:9)-deltaPos(deltaCt,4:6);
                
                % remember iter,success (better at beginning, but I don't want
                % to change all the indices now
                deltaPos(deltaCt,22:23) = debugData.trackResults(t,1,s).info(6:7);
                
                % remember sigma0
                if length(debugData.trackResults(t,1,s).sigma0)>1
                deltaPos(deltaCt,24) = debugData.trackResults(t,1,s).sigma0(2);
                end
            end

            % read idlisttrack if not estimated tag
            if idlisttrack(t).linklist(iTag,3) ~= 1
                trackCt = trackCt + 1;
                trackPos(trackCt,1) = t;
                trackPos(trackCt,2) = iPos;
                trackPos(trackCt,3) = length(idlisttrack(t).info.sourceList);

                trackPos(trackCt,4:6) = idlisttrack(t).linklist(iTag,[10,9,11])./pix2mu;
                trackPos(trackCt,10:12) = truePos - trackPos(trackCt,4:6);
                trackPos(trackCt,7:9) = truePos - floor(truePos);


                trackPos(trackCt,13:15) = trackQDiag( (iTag-1)*3+1:iTag*3 );

            end
        end
    end
end

% plot

% Figure 1: true error vs estimated error
% detector
figure('Name','accuracy vs. uncertainty')
sh = subplot(2,2,1);

ph=plot(sourcePos(:,9),sourcePos(:,12),'.r',...
    sourcePos(:,10),sourcePos(:,13),'.g',...
    sourcePos(:,11),sourcePos(:,14),'.b');
legend(ph,{'x','y','z'});
xlabel(sh,'accuracy')
ylabel(sh,'uncertainty')
title(sh,'detector')

% tracker(end)
sh = subplot(2,2,2);

ph=plot(trackPos(:,12),trackPos(:,15),'.b',...
    trackPos(:,10),trackPos(:,13),'.r',...
    trackPos(:,11),trackPos(:,14),'.g');
legend(ph,{'x','y','z'});
xlabel(sh,'accuracy')
ylabel(sh,'uncertainty')
title(sh,'tracker')

% delta
sh = subplot(2,2,3);

success = deltaPos(:,19)~=0;
ph=plot(deltaPos(success,16),deltaPos(success,19),'.r',...
    deltaPos(success,17),deltaPos(success,20),'.g',...
    deltaPos(success,18),deltaPos(success,21),'.b');
legend(ph,{'x','y','z'});
xlabel(sh,'accuracy')
ylabel(sh,'uncertainty')
title(sh,'delta')

% detector vs. tracker
sh = subplot(2,2,4);


ph=plot3(sourcePos(:,9),trackPos(:,10),sourcePos(:,9)-trackPos(:,10),'.r',...
    sourcePos(:,10),trackPos(:,11),sourcePos(:,10)-trackPos(:,11),'.g',...
    sourcePos(:,11),trackPos(:,12),sourcePos(:,11)-trackPos(:,12),'.b');
legend(ph,{'x','y','z'});
xlabel(sh,'detector')
ylabel(sh,'tracker')
zlabel(sh,'detector - tracker')
view(2)
title(sh,'detector vs. tracker')

% figure 2: sub-pixel position
% collect data
sourceIm = zeros(6,6,6);
trackIm = zeros(6,6,6);
deltaImS = zeros(6,6,6);
deltaImT = zeros(6,6,6);
[x, y, z] = ndgrid(1:6);
xs = (x(:)-1) * 0.1;
ys = (y(:)-1) * 0.1;
zs = (z(:)-1) * 0.1;
for i=1:216
    sppIdx = isApproxEqual(sourcePos(:,6),xs(i)) &...
        isApproxEqual(sourcePos(:,7), ys(i)) &...
        isApproxEqual(sourcePos(:,8) , zs(i));
    sourceIm(x(i),y(i),z(i)) = mean(sqrt(sum(sourcePos(sppIdx,9:11).^2,2)));
    
    sppIdx = isApproxEqual(deltaPos(:,10) , xs(i)) &...
        isApproxEqual(deltaPos(:,11) , ys(i)) &...
        isApproxEqual(deltaPos(:,12) , zs(i));
    deltaImS(x(i),y(i),z(i)) = mean(sqrt(sum(deltaPos(sppIdx,16:18).^2,2)));
    
    sppIdx = isApproxEqual(deltaPos(:,13) , xs(i)) &...
        isApproxEqual(deltaPos(:,14) , ys(i)) &...
        isApproxEqual(deltaPos(:,15) , zs(i));
    deltaImT(x(i),y(i),z(i)) = mean(sqrt(sum(deltaPos(sppIdx,16:18).^2,2)));
    
    sppIdx = isApproxEqual(trackPos(:,7) , xs(i)) &...
        isApproxEqual(trackPos(:,8) , ys(i)) &...
        isApproxEqual(trackPos(:,9) , zs(i));
    trackIm(x(i),y(i),z(i)) = mean(sqrt(sum(trackPos(sppIdx,10:12).^2,2)));
end

% figure('Name','uncertainty vs. sub-pixel-position')
% 
% % detector
% sh = subplot(2,2,1);
% 
% scatter3(xs(:), ys(:), zs(:), 1, sourceIm(:));
% xlabel(sh,'x')
% ylabel(sh,'y')
% zlabel(sh,'z')
% title(sh,'detector')

imarisShowArray(cat(5,sourceIm, deltaImS, deltaImT, trackIm));

% figure 3 delta xyz of tracker vs var vs nSources

figure('Name','accuracy vs uncertainty vs # of sources')
sh = subplot(2,2,1);
% x
ph=plot3(trackPos(:,10),trackPos(:,13),trackPos(:,3),'.r');
xlabel(sh,'accuracy')
ylabel(sh,'uncertainty')
zlabel(sh,'# of sources')
title(sh,'x')

sh = subplot(2,2,2);
% y
ph=plot3(trackPos(:,11),trackPos(:,14),trackPos(:,3),'.g');
xlabel(sh,'accuracy')
ylabel(sh,'uncertainty')
zlabel(sh,'# of sources')
title(sh,'y')


sh = subplot(2,2,3);
% z
ph=plot3(trackPos(:,12),trackPos(:,15),trackPos(:,3),'.b');
xlabel(sh,'accuracy')
ylabel(sh,'uncertainty')
zlabel(sh,'# of sources')
title(sh,'z')


% xyz
sh = subplot(2,2,4);

ph=plot3(sqrt(sum(trackPos(:,11:14).^2,2)),...
    sqrt(sum(trackPos(:,15).^2,2)),trackPos(:,3),'.k');
xlabel(sh,'accuracy')
ylabel(sh,'uncertainty')
zlabel(sh,'# of sources')
title(sh,'xyz')

% figure 4 delta delta vs. sigmaZero
figure('Name','accuracy (delta) vs. sigma0)')
plot(sqrt(sum(deltaPos(success,16:18).^2,2)),deltaPos(success,24),'.g',...
    sqrt(sum(deltaPos(~success,16:18).^2,2)),deltaPos(~success,24),'.g');

