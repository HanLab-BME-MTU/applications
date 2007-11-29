function cellDataAll = measureSpindles(idlist,dicMaxProj,doPlot)
%MEASURESPINDLES groups tags from several cells and measures spindle length
%
% SYNOPSIS: measureSpindles
%
% INPUT
%
% OUTPUT cellDataAll : nImages-by-3 cell array
%           1st row: for every cell in the image:
%               [nTags, spindleLength, smallest spb-cen-spb angle, 
%                   are spb-amps more similar to each other than to
%                   cen-amps?, mean(cenInt)/mean(spbInt),
%                   delta(cenInt)/mean(cenInt), delta(spbInt)/mean(spbInt)
%                   dist(spb-cen) x2, dist(cen1-cen2),
%                   normed projections of dist(spb-cen) x2]
%           2nd row: [#G1,#Pro,#M3,#M4]
%           3rd row: first two cols of 1st row for metaphase cells -> make
%                    into plotData for trapezoids
%
% REMARKS to add: cen-twist
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

def_doPlot = true;


if nargin == 0 || isempty(idlist)
    % load idlist
    idlistList = loadIdlistList([],[],true);
    
    % select files
selection = listSelectGUI({idlistList.name});

if isempty(selection)
    cellDataAll = [];
    return
end

% only keep selection
idlistList = idlistList(selection);
    
else
    idlistList.idlist = idlist;
    idlistList.dirName = pwd;
end



nIdlists = length(idlistList);

cellDataAll = cell(nIdlists,3);

% init doPlot
if nargin < 3 || isempty(doPlot)
doPlot = def_doPlot;
end

for iIdlist = 1:nIdlists

    idlistDir = idlistList(iIdlist).dirName;
    idlist = idlistList(iIdlist).idlist;

    if doPlot
        if (nargin < 2 || isempty(dicMaxProj))
            % load DIC image.

            % check for spindle image
            dicMaxProjName = dir(fullfile(idlistDir,'DIC_*.tif'));
            if ~isempty(dicMaxProjName)
                dicMaxProj=imread(fullfile(idlistDir,dicMaxProjName.name));
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
                    imwrite(dicMaxProj,[idlistDir,filesep,dicMaxProjName.name,'.tif'],'TIFF');
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
    % make clusters if the longest distance between two points is less than
    % 2 microns
    if doPlot
        figure,dendrogram(links,0,'Colorthreshold',2.2);
    end
    labels = cluster(links,'cutoff',2.2,'criterion','distance');

    % for every "cell": find longest distance, number of points, angle opposite
    % the longest side, flag whether the intensities of the two tags associated
    % with the longest side are more similar than the other tag(s)
    nCells = max(labels);
    cellData = zeros(nCells, 12);

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
        cenPosNorm = [0,0];
        


        if nTags > 2
            % find which tags the maxDist belongs to
            [row,col]=find(tril(ones(nTags),-1));
            % maxDistCidx: largest distance between tags - these are spbs
            maxDistCidx = cIdx([row(maxDistIdx);col(maxDistIdx)]);
            % minDistCidx: all the other tags
            minDistCidx = setdiff(cIdx,maxDistCidx);
            
            % maxAmp: spb amp (max b/c of maxDist)
            maxAmp = idlist(1).linklist(maxDistCidx,8);
            deltaMaxAmp = diff(maxAmp);
            % minAmp: cen-amp
            minAmp = idlist(1).linklist(minDistCidx,8);
            
            % if there are more than 4 tags, remove centromere with the
            % smallest amp
            [dummy,sortIdx] = sort(minAmp);
            killIdx = sortIdx(1:end-2);
            % also remove the index
            minDistCidx(killIdx) = [];
            minAmp(killIdx) = [];
            
            v_spb = coords(maxDistCidx(1),:) - coords(maxDistCidx(2),:);
            [n_spb,e_spb] = normList(v_spb);

            % for each minDistCidx: Find angle, similarAmplitude
            for i=length(minDistCidx):-1:1
                % angle - the angle spb1-cenX-spb2
                v(1,:) = coords(maxDistCidx(1),:)-coords(minDistCidx(i),:);
                v(2,:) = coords(maxDistCidx(2),:)-coords(minDistCidx(i),:);
                n = normList(v);
                angle(i) = 180/pi * acos(dot(v(1,:),v(2,:))/(n(1)*n(2)));

                % amplitude: spb amplitudes are considered similar, if
                % their difference is smaller than the difference between
                % any spindle pole and the centromere.
                similarAmplitude(i) =  deltaMaxAmp < max(abs(maxAmp-minAmp(i)));
                
                % distance to closest spb                
                [minDist2Spb(i),closeIdx] = min(n);
                
                % don't know what direction the vectors are, and don't care
                cenPosNorm(i) = abs(dot(v(closeIdx,:), e_spb, 2));
            end
            
            
            % cen relative to SPB, cenDelta rel to CEN
            relCenInt = mean(minAmp)/mean(maxAmp);
            if nTags > 3
                relCenDiff = diff(minAmp)/mean(minAmp);
                cenDist = distMat(coords(minDistCidx,:));
                cenDist = max(cenDist(:));
            else
                relCenDiff = 0;
            end
            relSpbDiff = diff(maxAmp)/mean(maxAmp);
            
            
            

        else
            if nTags == 1
                n_spb = 0;
            else
                % assume g1
                n_spb=1;
            end
            relCenInt = 0;
            relCenDiff = 0;
            relSpbDiff = 0;
        end

        cellData(c,:) = [nTags,n_spb,min(angle),min(similarAmplitude),relCenInt,relCenDiff,relSpbDiff,minDist2Spb(1:2),cenDist,cenPosNorm(1:2)];

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
    nM4 = sum(cellData(metaIdxL,1) > 3); % we only count 2 cens if nsp==5
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
        % plotData: spindleLength, bin# of spb-cen distance, weight,
        % movie#.cell#
        
        % collect lists
        plotData = zeros(nM3+2*nM4,4);
        metaIdx = find(metaIdxL);
        pdCt = 1;
        for idx=metaIdx';
            currentCell = cellData(idx,:);
            
            % find bin#:  d(s-c)/d(s-s) mapped onto 26 1/24 bins from
            % -1/24:0 to 24/24:25/24
            bin =ceil((currentCell(11:12)/currentCell(2)+1/24)*24);
            
            % remove outside bins
            if currentCell(1) == 3
                bin(2) = [];
            end
            bin(bin == 0) = 1;
            bin(bin < 1) = [];
            bin(bin > 26) = [];
            
            % count
            nBins = length(bin);
            
            % get movie/cell number - assume there will be below 1000 cells
            cellIdx = iIdlist + idx/1000; % movie.cell
            
            % write plotData
            endRange = pdCt + nBins - 1;
            plotData(pdCt:endRange,1) = currentCell(2); % n_spb
            plotData(pdCt:endRange,2) = bin(:); % bin
            plotData(pdCt:endRange,3) = 1/nBins; % weight
            plotData(pdCt:endRange,4) = cellIdx;
            
            % up counter
            pdCt = endRange + 1;
        end
        
        % remove superfluous entries
        plotData(pdCt:end,:) = [];
        
        %cellDataAll{iIdlist,3} = [cellData(metaIdxL,2),cellData(metaIdxL,1)];
        cellDataAll{iIdlist,3} = plotData;
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

% show 3vs 4 spots (sorry, it's another hack)
a=cat(1,cellDataAll{:,1});
a=a(a(:,1)>2 & a(:,2)>0.5,1:2);
range = 0.5:0.25:2;
c3=hist(a(a(:,1)==3,2),range);
c4=hist(a(a(:,1)==4,2),range);
figure('Name','3 (blue) and 4 (red) spots'),area(range,[c3',c4'])
xlabel('Spindle Length (\mum)')
figure('Name','3 (blue) and 4 (red) spots'),area(range,[c3',c4']./repmat(nansum([c3',c4'],2),1,2))
xlabel('Spindle Length (\mum)')
% reset warnings
warning(warningState);

% Eugenio's request of Oct 16 2007: plot relCenInt for 3 and 4 spots
figure('Name','relative centromere intensity (3=b/4=r)')
plot(data1(d31,2),data1(d31,5),'.b',data1(d41,2),data1(d41,5),'.r')
xlabel('Spindle Length (\mum)')
ylabel('$$\frac{\Sigma(cenInt)/numCen}{\Sigma(spbInt)/2}$$','Interpreter','latex')