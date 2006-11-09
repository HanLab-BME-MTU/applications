function cellDataAll = measureSpindles(idlist,dicMaxProj)
%MEASURESPINDLES groups tags from several cells and measures spindle length
%
% SYNOPSIS: measureSpindles
%
% INPUT
%
% OUTPUT
%
% REMARKS
%
% created with MATLAB ver.: 7.2.0.232 (R2006a) on Windows_NT
%
% created by: Jonas Dorn
% DATE: 28-Sep-2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% remember warning state
warningState = warning;
warning off Images:initSize:adjustingMag
warning off MATLAB:intConvertNonIntVal
warning off MATLAB:divideByZero


if nargin == 0 || isempty(idlist)
    % load idlist
    idlistList = loadIdlistList([],[],true);
else
    idlistList.idlist = idlist;
    idlistList.dirName = pwd;
end

nIdlists = length(idlistList);

cellDataAll = cell(nIdlists,3);

% init doPlot
doPlot = true;

for iIdlist = 1:nIdlists

    idlistDir = idlistList(iIdlist).dirName;
    idlist = idlistList(iIdlist).idlist;

    if doPlot
        if (nargin < 2 || isempty(dicMaxProj))
            % load DIC image.

            % check for spindle image
            dicMaxProjName = dir(fullfile(idlistDir,'DIC_*.tif'));
            if ~isempty(dicMaxProjName)
                dicMaxProj=imread(dicMaxProjName.name);
                doPlot = 1;
            else
                % try to find DIC
                dicMaxProjName = dir(fullfile(idlistDir,'DIC*'));
                if isempty(dicMaxProjName)
                    warning('no DIC image found')
                    doPlot = 0;
                else
                    % load DIC image
                    dicImage=r3dread(fullfile(idlistDir,dicMaxProjName(1).name));
                    % max projection
                    dicMaxProj = max(dicImage,[],3);
                    % save maxprojection
                    imwrite(dicMaxProj,[dicMaxProjName.name,'.tif'],'TIFF');
                    clear dicImage
                    doPlot = 1;
                end
            end

        elseif length(dicMaxProj) == 1
            doPlot = dicMaxProj;
        else
            doPlot = 1;
        end
    end

    % calculate distances
    coords = idlist(1).linklist(:,9:11);
    distList = pdist(coords); % default is euclidean

    % cluster with complete linkage
    links = linkage(distList,'complete');
    % make clusters if the longest distance between two points is less than 3
    % microns
    if doPlot
        figure,dendrogram(links,0,'Colorthreshold',2.5);
    end
    labels = cluster(links,'cutoff',2.5,'criterion','distance');

    % for every "cell": find longest distance, number of points, angle opposite
    % the longest side, flag whether the intensities of the two tags associated
    % with the longest side are more similar than the other tag(s)
    nCells = max(labels);
    cellData = zeros(nCells, 10);

    for c = 1:nCells
        % init vars
        angle = 180;
        similarAmplitude = 1;


        cIdx = find(labels == c);
        nTags = length(cIdx);
        dist = pdist(coords(cIdx,:));
        [maxDist,maxDistIdx] = max(dist);

        minDist2Spb = [0,0];
        cenDist = 0;
        


        if nTags > 2
            % find which tags the maxDist belongs to
            [row,col]=find(tril(ones(nTags),-1));
            maxDistCidx = cIdx([row(maxDistIdx);col(maxDistIdx)]);
            minDistCidx = setdiff(cIdx,maxDistCidx);

            maxAmp = idlist(1).linklist(maxDistCidx,8);
            deltaMaxAmp = diff(maxAmp);
            minAmp = idlist(1).linklist(minDistCidx,8);
            

            % for each minDistCidx: Find angle, similarAmplitude
            for i=length(minDistCidx):-1:1
                % angle
                v1 = coords(maxDistCidx(1),:)-coords(minDistCidx(i),:);
                v2 = coords(maxDistCidx(2),:)-coords(minDistCidx(i),:);
                angle(i) = 180/pi * acos(dot(v1,v2)/(norm(v1)*norm(v2)));

                % amplitude
                similarAmplitude(i) =  deltaMaxAmp < max(maxAmp-minAmp(i));
                
                % distance to closest spb
                minDist2Spb(i) = min(distMat2(coords(maxDistCidx,:),coords(minDistCidx(i),:)));
            end

            % cen relative to SPB, cenDelta rel to CEN
            relCenInt = mean(minAmp)/mean(maxAmp);
            if nTags == 4
                relCenDiff = diff(minAmp)/mean(minAmp);
                cenDist = distMat(coords(minDistCidx,:));
                cenDist = max(cenDist(:));
            else
                relCenDiff = 0;
            end
            relSpbDiff = diff(maxAmp)/mean(maxAmp);

        else
            if nTags == 1
                dist = 0;
            end
            relCenInt = 0;
            relCenDiff = 0;
            relSpbDiff = 0;
        end

        cellData(c,:) = [nTags,max(dist),min(angle),min(similarAmplitude),relCenInt,relCenDiff,relSpbDiff,minDist2Spb(1:2),cenDist];

    end

    if doPlot
        figure('Name',idlist(1).stats.name),imshow(dicMaxProj',[]);
        hold on
        for c = 1:nCells
            cIdx = find(labels == c);
            xy = coords(cIdx,[2,1]);
            if isEven(c)
                plot(xy(:,1)/0.0663,xy(:,2)/0.066,'o','MarkerSize',4,'Color',extendedColors(c/2));
            else
                plot(xy(:,1)/0.0663,xy(:,2)/0.066,'s','MarkerSize',4,'Color',extendedColors((c+1)/2));
            end
        end
    end

    % display statistics
    % metaphase: minimum of 1um spindle length, similar spb-amp
    metaIdxL = cellData(:,1) == 4 | (cellData(:,1) == 3 & cellData(:,4) == 1 & cellData(:,2) > 1);
    proMetaIdxL = cellData(:,1) == 3 & ~metaIdxL;
    g1IdxL = cellData(:,1) <= 2;
    nM3 = sum(cellData(metaIdxL,1) == 3);
    nM4 = sum(cellData(metaIdxL,1) == 4);
    nPro = sum(proMetaIdxL);
    nG1 = sum(g1IdxL);
    nM = sum(metaIdxL);

    % if any(metaIdxL)
    %     figure('Name',idlist.name)
    %     boxplot(cellData(metaIdxL,2),cellData(metaIdxL,1),'notch','on')
    % end

    disp(sprintf('%1.2f%% G1, %1.2f%% proM, %1.2f%% M, %1.2f%% 3spM, %1.2f%% 4spM, Ntot %i, NM %i',...
        nG1/nCells*100,nPro/nCells*100,nM/nCells*100,nM3/nM*100,nM4/nM*100,nCells,nM))

    % collect statistics
    cellDataAll{iIdlist,1} = cellData;
    cellDataAll{iIdlist,2} = [nG1,nPro,nM3,nM4];
    if any(metaIdxL)
        cellDataAll{iIdlist,3} = [cellData(metaIdxL,2),cellData(metaIdxL,1)];
    end

    % reset DIC-image
    dicMaxProj = [];
end


data1=cat(1,cellDataAll{:,1});
metaIdxL = data1(:,1) == 4 | (data1(:,1) == 3 & data1(:,4) == 1 & data1(:,2) > 1);
    proMetaIdxL = data1(:,1) == 3 & ~metaIdxL;
    g1IdxL = data1(:,1) <= 2;
    nM3 = sum(data1(metaIdxL,1) == 3);
    nM4 = sum(data1(metaIdxL,1) == 4);
    nPro = sum(proMetaIdxL);
    nG1 = sum(g1IdxL);
    nM = sum(metaIdxL);
    nCells = size(data1,1);

    % if any(metaIdxL)
    %     figure('Name',idlist.name)
    %     boxplot(cellData(metaIdxL,2),cellData(metaIdxL,1),'notch','on')
    % end

    disp(sprintf('%1.2f%% G1, %1.2f%% proM, %1.2f%% M, %1.2f%% 3spM, %1.2f%% 4spM, Ntot %i, NM %i',...
        nG1/nCells*100,nPro/nCells*100,nM/nCells*100,nM3/nM*100,nM4/nM*100,nCells,nM))



% plot overall


% intStats for 4 spots
d41=data1(:,1)==4;
figure('Name','IntStats 4sp'),
[H,AX,BigAx,P] = plotmatrix(abs(data1(d41,[2,5:7,10])),'.');
set(get(AX(1,1),'yLabel'),'String','SpbDist (\mum)')
set(get(AX(2,1),'yLabel'),'String','relCenInt')
set(get(AX(3,1),'yLabel'),'String','relCenDiff')
set(get(AX(4,1),'yLabel'),'String','relSpbDiff')
set(get(AX(5,1),'yLabel'),'String','cenDist')
set(get(AX(5,1),'xLabel'),'String','SpbDist (\mum)')
set(get(AX(5,2),'xLabel'),'String','relCenInt')
set(get(AX(5,3),'xLabel'),'String','relCenDiff')
set(get(AX(5,4),'xLabel'),'String','relSpbDiff')
set(get(AX(5,5),'xLabel'),'String','cenDist')

figure,plot3([data1(d41,2);data1(d41,2)],abs([data1(d41,6);data1(d41,6)]),[data1(d41,10);data1(d41,10)],'.')

%figure,plot3([data1(d41,2);data1(d41,2)],abs([data1(d41,6);data1(d41,6)]),[data1(d41,8);data1(d41,9)],'.')

% intStats for 3 spots
d31=(data1(:,1) == 3 & data1(:,4) == 1 & data1(:,2) > 1);
figure('Name','IntStats 3sp'),
[H,AX,BigAx,P] = plotmatrix(abs(data1(d31,[2,5,7])),'.');
set(get(AX(1,1),'yLabel'),'String','SpbDist (\mum)')
set(get(AX(2,1),'yLabel'),'String','relCenInt')
set(get(AX(3,1),'yLabel'),'String','relSpbDiff')
set(get(AX(3,1),'xLabel'),'String','SpbDist (\mum)')
set(get(AX(3,2),'xLabel'),'String','relCenInt')
set(get(AX(3,3),'xLabel'),'String','relSpbDiff')




% reset warnings
warning(warningState);
