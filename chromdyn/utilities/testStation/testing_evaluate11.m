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

if nargin == 0
    % load from file
    isLaptop = ispc && strcmp(getenv('HOME'),'G:\matlab');

    if isLaptop
    dataPath = 'G:\matlab\trackTest';
else
    [bioDataDir, oldDir] = cdBiodata(0);
    % go to right dir
    cd(['..' filesep 'trackTest']);
    dataPath = pwd;
    cd(oldDir);
    end
    % cd to directory. it's easier
cd(sprintf('%s\\testMovie_902',dataPath))

[dummy, tokens] = searchFiles(...
    'debugData_P902_S([\d.]+)_G([\d]+)_i([\d]+)');
% make a matrix of numbers
sgiIdx = str2double(tokens);

% sort
[sgiIdx,sortIdx] = sortrows(sgiIdx,[-1,2,3]);

% find number of different snrs, gradients, repeats
nSnr = length(unique(sgiIdx(:,1)));
nGrad = length(unique(sgiIdx(:,2)));
nRepeats = length(unique(sgiIdx(:,3)));

% loop to load idlisttrack and debugData
rawData(1:nSnr,1:nGrad,1:nRepeats) = struct('debugData',[],'idlisttrack',[]);
allData(1:nSnr,1:nGrad,1:nRepeats) = struct('snr',[],'grad',[]);

for ct = 1:nSnr*nGrad*nRepeats
    % load debugData, idlisttrack. remove fstats for now
            load(sprintf('debugData_P902_S%1.2f_G%i_i%i.mat',sgiIdx(ct,:)))
            load(sprintf('idlisttrack_NS5_S%1.2f_G%i_i%i.mat',sgiIdx(ct,:)))
            debugData.fStats = [];
            rawData(ct).debugData = debugData;
            rawData(ct).idlisttrack = idlisttrack;
            allData(ct).snr = sgiIdx(ct,1);
            allData(ct).grad = sgiIdx(ct,2);
end

% load dataProperties
load dataProperties



else
    rawData.debugData = debugData;
    rawData.idlisttrack = idlisttrack;

    allData.snr = -1;
    allData.grad = -1;

    nSnr = 1;
    nGrad = 1;
end

% remove debugData from memory
clear debugData

% prepare collecting data
positions = loadPositions(902);
positions = positions(:,[1,2,3],:);
nTags = size(positions,3);

% loop data
for iSnr = 1:nSnr
    for iGrad = 1:nGrad



        pix2mu = [dataProperties.PIXELSIZE_XY,dataProperties.PIXELSIZE_XY,dataProperties.PIXELSIZE_Z];



        id2posIdx = zeros(nTags,1);
        for i=1:nTags
            id2posIdx(i) = str2double(rawData(iSnr,iGrad).idlisttrack(1).stats.labelcolor{i});
        end


        numSources = size(rawData(iSnr,iGrad).debugData.trackResults,3);
        nTimepoints = length(rawData(iSnr,iGrad).idlisttrack);

        % sourcePos:
        % t, iTag, source XYZ, subpixelposition XYZ (groundTruth), sourceDelta XYZ,
        % sqrt(sourceVar XYZ)

        % deltaPos
        % t, iTag, tSource, delta, trueDelta, sourceSubPixel, targetSubPixel, deltaDelta,
        % sqrt(deltaVar), iter, success, sigma0_end

        % trackPos
        % t, iTag, nSources, pos, true pos (subpixel), deltaPos, sqrt(var)

        allData(iSnr,iGrad).sourcePos = zeros(nTimepoints * nTags,14);
        allData(iSnr,iGrad).deltaPos = zeros(nTimepoints * nTags * numSources,24);
        allData(iSnr,iGrad).trackPos = zeros(nTimepoints * nTags,15);
        sourceCt = 0;
        deltaCt = 0;
        trackCt = 0;


        for t=1:nTimepoints
            if ~isempty(rawData(iSnr,iGrad).idlisttrack(t).linklist)
                trackQDiag = diag(rawData(iSnr,iGrad).idlisttrack(t).info.totalQ_Pix);
                for iTag=1:nTags
                    iPos = id2posIdx(iTag);

                    truePos = positions(t,:,iPos);


                    if isempty(rawData(iSnr,iGrad).debugData.trackResults(t,1,1).info) ||...
                            rawData(iSnr,iGrad).debugData.trackResults(t,1,1).info(2) == 1;
                        % it's a source
                        sourceCt = sourceCt + 1;

                        allData(iSnr,iGrad).sourcePos(sourceCt,1) = t;
                        allData(iSnr,iGrad).sourcePos(sourceCt,2) = iPos;
                        allData(iSnr,iGrad).sourcePos(sourceCt,3:5) = rawData(iSnr,iGrad).debugData.trackResults(t,iTag,1).sourcePos;
                        allData(iSnr,iGrad).sourcePos(sourceCt,12:14) = sqrt(rawData(iSnr,iGrad).debugData.trackResults(t,iTag,1).sourceVar);

                        allData(iSnr,iGrad).sourcePos(sourceCt,9:11) = truePos - allData(iSnr,iGrad).sourcePos(sourceCt,3:5);
                        allData(iSnr,iGrad).sourcePos(sourceCt,6:8) = truePos - floor(truePos);

                        % start looking at start/end from #2
                        s1=2;
                    else
                        s1=1;
                    end

                    % loop tracking
                    for s = s1:numSources
                        deltaCt = deltaCt + 1;
                        allData(iSnr,iGrad).deltaPos(deltaCt,1) = t;
                        allData(iSnr,iGrad).deltaPos(deltaCt,2) = iPos;

                        % if clause for error correction
                        if size(rawData(iSnr,iGrad).debugData.trackResults(t,iTag,s).startEndDelta,2) == 6
                            allData(iSnr,iGrad).deltaPos(deltaCt,4:6) = rawData(iSnr,iGrad).debugData.trackResults(t,iTag,s).startEndDelta(2,4:6);
                        elseif size(rawData(iSnr,iGrad).debugData.trackResults(t,iTag,s).startEndDelta,1) == 6
                            allData(iSnr,iGrad).deltaPos(deltaCt,4:6) = rawData(iSnr,iGrad).debugData.trackResults(t,iTag,s).startEndDelta(4:6);
                        else
                            allData(iSnr,iGrad).deltaPos(deltaCt,4:6) = rawData(iSnr,iGrad).debugData.trackResults(t,iTag,s).startEndDelta(2,:);
                        end
                        allData(iSnr,iGrad).deltaPos(deltaCt,19:21) = sqrt(rawData(iSnr,iGrad).debugData.trackResults(t,iTag,s).deltaVar);

                        sourceT = rawData(iSnr,iGrad).debugData.trackResults(t,1,s).info(3);
                        allData(iSnr,iGrad).deltaPos(deltaCt,3) = sourceT;

                        truePosS = positions(sourceT,:,iPos);

                        allData(iSnr,iGrad).deltaPos(deltaCt,7:9) = truePos - truePosS;
                        allData(iSnr,iGrad).deltaPos(deltaCt,10:12) = truePosS - floor(truePosS);
                        allData(iSnr,iGrad).deltaPos(deltaCt,13:15) = truePos - floor(truePos);
                        allData(iSnr,iGrad).deltaPos(deltaCt,16:18) = allData(iSnr,iGrad).deltaPos(deltaCt,7:9)-allData(iSnr,iGrad).deltaPos(deltaCt,4:6);

                        % remember iter,success (better at beginning, but I don't want
                        % to change all the indices now
                        allData(iSnr,iGrad).deltaPos(deltaCt,22:23) = rawData(iSnr,iGrad).debugData.trackResults(t,1,s).info(6:7);

                        % remember sigma0
                        if length(rawData(iSnr,iGrad).debugData.trackResults(t,1,s).sigma0)>1
                            allData(iSnr,iGrad).deltaPos(deltaCt,24) = rawData(iSnr,iGrad).debugData.trackResults(t,1,s).sigma0(2);
                        end
                    end

                    % read idlisttrack if not estimated tag
                    if rawData(iSnr,iGrad).idlisttrack(t).linklist(iTag,3) ~= 1
                        trackCt = trackCt + 1;
                        allData(iSnr,iGrad).trackPos(trackCt,1) = t;
                        allData(iSnr,iGrad).trackPos(trackCt,2) = iPos;
                        allData(iSnr,iGrad).trackPos(trackCt,3) = length(rawData(iSnr,iGrad).idlisttrack(t).info.sourceList);

                        allData(iSnr,iGrad).trackPos(trackCt,4:6) = rawData(iSnr,iGrad).idlisttrack(t).linklist(iTag,[10,9,11])./pix2mu;
                        allData(iSnr,iGrad).trackPos(trackCt,10:12) = truePos - allData(iSnr,iGrad).trackPos(trackCt,4:6);
                        allData(iSnr,iGrad).trackPos(trackCt,7:9) = truePos - floor(truePos);


                        allData(iSnr,iGrad).trackPos(trackCt,13:15) = trackQDiag( (iTag-1)*3+1:iTag*3 );

                    end
                end
            end
        end

        % plot

        % Figure 1: true error vs estimated error
        % detector
        figure('Name',sprintf('S %1.2f, G %i, accuracy vs. uncertainty',allData(iSnr,iGrad).snr, allData(iSnr,iGrad).grad))
        sh = subplot(2,2,1);

        ph=plot(allData(iSnr,iGrad).sourcePos(:,9),allData(iSnr,iGrad).sourcePos(:,12),'.r',...
            allData(iSnr,iGrad).sourcePos(:,10),allData(iSnr,iGrad).sourcePos(:,13),'.g',...
            allData(iSnr,iGrad).sourcePos(:,11),allData(iSnr,iGrad).sourcePos(:,14),'.b');
        legend(ph,{'x','y','z'});
        xlabel(sh,'accuracy')
        ylabel(sh,'uncertainty')
        title(sh,'detector')

        % tracker(end)
        sh = subplot(2,2,2);

        ph=plot(allData(iSnr,iGrad).trackPos(:,12),allData(iSnr,iGrad).trackPos(:,15),'.b',...
            allData(iSnr,iGrad).trackPos(:,10),allData(iSnr,iGrad).trackPos(:,13),'.r',...
            allData(iSnr,iGrad).trackPos(:,11),allData(iSnr,iGrad).trackPos(:,14),'.g');
        legend(ph,{'x','y','z'});
        xlabel(sh,'accuracy')
        ylabel(sh,'uncertainty')
        title(sh,'tracker')

        % delta
        sh = subplot(2,2,3);

        success = allData(iSnr,iGrad).deltaPos(:,19)~=0;
        ph=plot(allData(iSnr,iGrad).deltaPos(success,16),allData(iSnr,iGrad).deltaPos(success,19),'.r',...
            allData(iSnr,iGrad).deltaPos(success,17),allData(iSnr,iGrad).deltaPos(success,20),'.g',...
            allData(iSnr,iGrad).deltaPos(success,18),allData(iSnr,iGrad).deltaPos(success,21),'.b');
        legend(ph,{'x','y','z'});
        xlabel(sh,'accuracy')
        ylabel(sh,'uncertainty')
        title(sh,'delta')

        % detector vs. tracker
        sh = subplot(2,2,4);


        ph=plot3(allData(iSnr,iGrad).sourcePos(:,9),allData(iSnr,iGrad).trackPos(:,10),allData(iSnr,iGrad).sourcePos(:,9)-allData(iSnr,iGrad).trackPos(:,10),'.r',...
            allData(iSnr,iGrad).sourcePos(:,10),allData(iSnr,iGrad).trackPos(:,11),allData(iSnr,iGrad).sourcePos(:,10)-allData(iSnr,iGrad).trackPos(:,11),'.g',...
            allData(iSnr,iGrad).sourcePos(:,11),allData(iSnr,iGrad).trackPos(:,12),allData(iSnr,iGrad).sourcePos(:,11)-allData(iSnr,iGrad).trackPos(:,12),'.b');
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
            sppIdx = isApproxEqual(allData(iSnr,iGrad).sourcePos(:,6),xs(i)) &...
                isApproxEqual(allData(iSnr,iGrad).sourcePos(:,7), ys(i)) &...
                isApproxEqual(allData(iSnr,iGrad).sourcePos(:,8) , zs(i));
            sourceIm(x(i),y(i),z(i)) = mean(sqrt(sum(allData(iSnr,iGrad).sourcePos(sppIdx,9:11).^2,2)));

            sppIdx = isApproxEqual(allData(iSnr,iGrad).deltaPos(:,10) , xs(i)) &...
                isApproxEqual(allData(iSnr,iGrad).deltaPos(:,11) , ys(i)) &...
                isApproxEqual(allData(iSnr,iGrad).deltaPos(:,12) , zs(i));
            deltaImS(x(i),y(i),z(i)) = mean(sqrt(sum(allData(iSnr,iGrad).deltaPos(sppIdx,16:18).^2,2)));

            sppIdx = isApproxEqual(allData(iSnr,iGrad).deltaPos(:,13) , xs(i)) &...
                isApproxEqual(allData(iSnr,iGrad).deltaPos(:,14) , ys(i)) &...
                isApproxEqual(allData(iSnr,iGrad).deltaPos(:,15) , zs(i));
            deltaImT(x(i),y(i),z(i)) = mean(sqrt(sum(allData(iSnr,iGrad).deltaPos(sppIdx,16:18).^2,2)));

            sppIdx = isApproxEqual(allData(iSnr,iGrad).trackPos(:,7) , xs(i)) &...
                isApproxEqual(allData(iSnr,iGrad).trackPos(:,8) , ys(i)) &...
                isApproxEqual(allData(iSnr,iGrad).trackPos(:,9) , zs(i));
            trackIm(x(i),y(i),z(i)) = mean(sqrt(sum(allData(iSnr,iGrad).trackPos(sppIdx,10:12).^2,2)));
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

        % imarisShowArray(cat(5,sourceIm, deltaImS, deltaImT, trackIm));

        % figure 3 delta xyz of tracker vs var vs nSources

        figure('Name',sprintf('S %1.2f, G %i, accuracy vs uncertainty vs # of sources',allData(iSnr,iGrad).snr, allData(iSnr,iGrad).grad))
        sh = subplot(2,2,1);
        % x
        ph=plot3(allData(iSnr,iGrad).trackPos(:,10),allData(iSnr,iGrad).trackPos(:,13),allData(iSnr,iGrad).trackPos(:,3),'.r');
        xlabel(sh,'accuracy')
        ylabel(sh,'uncertainty')
        zlabel(sh,'# of sources')
        title(sh,'x')

        sh = subplot(2,2,2);
        % y
        ph=plot3(allData(iSnr,iGrad).trackPos(:,11),allData(iSnr,iGrad).trackPos(:,14),allData(iSnr,iGrad).trackPos(:,3),'.g');
        xlabel(sh,'accuracy')
        ylabel(sh,'uncertainty')
        zlabel(sh,'# of sources')
        title(sh,'y')


        sh = subplot(2,2,3);
        % z
        ph=plot3(allData(iSnr,iGrad).trackPos(:,12),allData(iSnr,iGrad).trackPos(:,15),allData(iSnr,iGrad).trackPos(:,3),'.b');
        xlabel(sh,'accuracy')
        ylabel(sh,'uncertainty')
        zlabel(sh,'# of sources')
        title(sh,'z')


        % xyz
        sh = subplot(2,2,4);

        ph=plot3(sqrt(sum(allData(iSnr,iGrad).trackPos(:,11:14).^2,2)),...
            sqrt(sum(allData(iSnr,iGrad).trackPos(:,15).^2,2)),allData(iSnr,iGrad).trackPos(:,3),'.k');
        xlabel(sh,'accuracy')
        ylabel(sh,'uncertainty')
        zlabel(sh,'# of sources')
        title(sh,'xyz')

        % figure 4 delta delta vs. sigmaZero
        figure('Name',sprintf('S %1.2f, G %i, accuracy (delta) vs. sigma0)',allData(iSnr,iGrad).snr, allData(iSnr,iGrad).grad))
        plot(sqrt(sum(allData(iSnr,iGrad).deltaPos(success,16:18).^2,2)),allData(iSnr,iGrad).deltaPos(success,24),'.g',...
            sqrt(sum(allData(iSnr,iGrad).deltaPos(~success,16:18).^2,2)),allData(iSnr,iGrad).deltaPos(~success,24),'.g');

        
        %==========
        
        % remove data from memory
        rawData(iSnr,iGrad).debugData = [];
        rawData(iSnr,iGrad).idlisttrack = [];
    end % loop iGrad
end % loop iSNR

figure,
