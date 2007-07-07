function dataStruct = makiFitPlane(dataStruct,verbose)
%MAKIFITPLANE attempts to fit planes into kinetochore clusters
%
% SYNOPSIS: makiFitPlane(dataStruct)
%
% INPUT dataStruct (opt): maki data structure. If empty, it will be loaded
%                         via GUI
%       verbose (opt)   : 1: plot results in the end (default)
%                         2: also plot plane fit for every frame
%
% OUTPUT dataStruct.planeFit structure of length nTimepoints. Except for
%                   eigenVectors and eigenValues, fields are only filled
%                   for metaphase timepoints
%               .plane [a b c d] for the plane equation ax+by+cz-d=0
%               .planeCoord      spot coordinates relative to the plane.
%                   The first dimension is perpendicular to the plane, the
%                   second direction is the intersect between the plane and
%                   the xy-plane
%               .planeVectors    unit vectors of the plane coordinate
%                   system
%               .eigenVectors    eigenVectors of the spot covariance matrix
%               .eigenValues     corresponding eigenValues
%               .inlierIdx       index into initCoord of all spots that are
%                   considered to be inliers
%               .laggingIdx      index into initCoord of all spots that are
%                   considered to belong to lagging chromosomes
%               .metaphaseFrames frames where the cell is considered to be
%                   in metaphase
%               .distParms       [variance,skewness,kurtosis,pNormal]' for
%                   the planeCoord. pNormal is the p-value that the data
%                   comes from a normal distribution. The tabulated values
%                   for the Lilliefors test only go from 0.01 to 0.5!
%               .deltaP          p-value of the ks-test comparing the
%                   distribution of the distances of inlier spots from the
%                   metaphase plate of t and t+1
%               .deltaAngle      difference in angle between planes in t
%                   and t+1
%
%
% REMARKS
%
% created with MATLAB ver.: 7.4.0.287 (R2007a) on Windows_NT
%
% created by: jdorn
% DATE: 03-Jul-2007
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warningState = warning;
warning off stats:lillietest:OutOfRangeP

%TEST input
if nargin == 0 || isempty(dataStruct)
    dataStruct = makiLoadDataFile;
end
if nargin < 2 || isempty(verbose)
    verbose = 1;
end

% get coordinates, dataProperties, etc
initCoord = dataStruct.initCoord;
nTimepoints = length(initCoord);
nSpots = cat(1,initCoord.nSpots);

planeFit(1:nTimepoints) = struct('plane',[],'planeCoord',[],...
    'planeVectors',[],'eigenVectors',[],'eigenValues',[],...
    'inlierIdx',[],'laggingIdx',[],'metaphaseFrames',[],...
    'distParms',[],'deltaP',[],'deltaAngle',[]);

% loop through timepoints. Get covariance of point cloud, and the
% corresponding eigenvalues. Compare the two most similar eigenvalues to
% the third to find out how much of a plane there is. Consider a ratio of
% below 1 to roughly be metaphase
eigenValues = zeros(nTimepoints, 3);
eigenVectors = zeros(3,3,nTimepoints);  %vectors in cols
meanCoord = zeros(nTimepoints,3);

for t=1:nTimepoints
    % get data for eigenRatio, subsequent plane fitting
    [eigenVectors(:,:,t), eigenValues(t,:), meanCoord(t,:)] = ...
        eigenCalc(initCoord(t).allCoord(:,1:3));
    planeFit(t).eigenVectors = eigenVectors(:,:,t);
    planeFit(t).eigenValues = eigenValues(t,:);
end

% metaphaseFrames are those with eigenRatio<1
eigenRatio = eigenValues(:,1)./mean(eigenValues(:,2:3),2);
goodFramesL = (eigenRatio(:,1)<1);
% often, the eigenRatio starts low and then peaks before going back down
% find this by using bwlabel and identify the largest group
goodFramesLb = bwlabel(goodFramesL);
[idx,cts] = countEntries(goodFramesLb);
goodLabel = idx(find(cts(2:end)==max(cts(2:end)))+1);
metaphaseFrames = find(goodFramesLb==goodLabel);

% loop metaphaseFrames to fit plane. Basically, the plane is defined by the
% eigenvectors and the average coordinate. However, there could be
% outliers. To identify outliers and improve the plane fit, the following
% strategy is chosen:
% - calculate distances from the plane
% - identify potential outliers via robustMean
% - refit the plane, using only inliers
% - identify potential outliers again (considering all points)
% - if no change in the outliers, exit the loop. Otherwise, continue


% store metaphaseFrames in pf1
planeFit(1).metaphaseFrames = metaphaseFrames;

% do only for good frames
for t = metaphaseFrames'
    done = false;
    % initially: assume no bad spots. Allow for 10 iterations
    badSpotIdxLOld = false(nSpots(t),10);
    ct = 1;

    while ~done

        % get distance from plane, in-plane coordinates by transformation
        % with inverse of in-plane vectors
        normal = eigenVectors(:,1,t);
        e_plane = zeros(3);
        e_plane(:,1) = normal;
        e_plane(:,2) = [-normal(2),normal(1),0]./sqrt(sum(normal(1:2).^2));
        e_plane(:,3) = cross(e_plane(:,1),e_plane(:,2));
        % planeCoord: [d,xplane,yplane]
        planeFit(t).planeCoord = ...
            (inv(e_plane)*...
            (initCoord(t).allCoord(:,1:3)-repmat(meanCoord(t,:),nSpots(t),1))')';

        % find outliers
        [dummy, dummy, goodSpotIdx] = robustMean(planeFit(t).planeCoord(:,1));
        badSpotIdxL = true(initCoord(t).nSpots,1);
        badSpotIdxL(goodSpotIdx) = false;

        %         % DEBUG plot plane in matlab
        %         pc=planeFit(t).planeCoord;
        %         pos = pc(:,1)>0;
        %         neg = pc(:,1)<0;
        %         figure,plot3(pc(pos,1),pc(pos,2),pc(pos,3),'.k',...
        %             pc(neg,1),pc(neg,2),pc(neg,3),'or',...
        %             pc(badSpotIdxL,1),pc(badSpotIdxL,2),pc(badSpotIdxL,3),'+b')

        % check whether there was any change
        if any(all(repmat(badSpotIdxL,1,10) == badSpotIdxLOld,1)) || ct == 10
            % done. Fill information into planeFit-structure

            planeFit(t).plane = [normal',meanCoord(t,:)*normal];
            planeFit(t).planeVectors = e_plane;
            planeFit(t).eigenVectors = eigenVectors(:,:,t);
            planeFit(t).eigenValues = eigenValues(t,:);

            % lagging chromosomes are outliers (until we can identify pairs
            planeFit(t).laggingIdx = find(badSpotIdxL);
            planeFit(t).inlierIdx = goodSpotIdx;

            % distribution parameters (do for all unit directions -
            % the second vector is also interesting, as it lies in the xy
            % plane, in which the metaphase plate should not be cut off
            % distribution parameters (rows):
            % var
            % skew
            % kurtosis
            % p for normal distribution (lilliefors test)
            % correct all parameters for bias
            planeFit(t).distParms = zeros(4,3);
            planeFit(t).distParms(1,:) = var(planeFit(t).planeCoord(goodSpotIdx,:));
            planeFit(t).distParms(2,:) = skewness(planeFit(t).planeCoord(goodSpotIdx,:),0);
            planeFit(t).distParms(3,:) = kurtosis(planeFit(t).planeCoord(goodSpotIdx,:),0);
            [dummy,planeFit(t).distParms(4,1)] = ...
                lillietest(planeFit(t).planeCoord(goodSpotIdx,1));
            [dummy,planeFit(t).distParms(4,2)] = ...
                lillietest(planeFit(t).planeCoord(goodSpotIdx,2));
            [dummy,planeFit(t).distParms(4,3)] = ...
                lillietest(planeFit(t).planeCoord(goodSpotIdx,3));

            % plot plane in matlab
            if verbose > 1
                [ygrid,zgrid] = meshgrid(...
                    linspace(min(planeFit(t).planeCoord(:,2)),...
                    max(planeFit(t).planeCoord(:,2)),5), ...
                    linspace(min(planeFit(t).planeCoord(:,3)),...
                    max(planeFit(t).planeCoord(:,3)),5));
                xgrid = zeros(5,5);
                pc=planeFit(t).planeCoord;
                pos = pc(:,1)>0;
                neg = pc(:,1)<0;
                figure('Name',...
                    sprintf('Metaphase plate frame %i for %s',...
                    t,dataStruct.projectName))
                plot3(pc(pos,1),pc(pos,2),pc(pos,3),'.k',...
                    pc(neg,1),pc(neg,2),pc(neg,3),'or',...
                    pc(badSpotIdxL,1),pc(badSpotIdxL,2),pc(badSpotIdxL,3),'+b')
                hold on
                mesh(xgrid,ygrid,zgrid,'EdgeColor',[0 0 1],'FaceAlpha',0);
                grid on
            end


            done = true;

        else
            % re-"fit" the plane. Update eigenVectors etc.
            [eigenVectors(:,:,t), eigenValues(t,:), meanCoord(t,:)] = ...
                eigenCalc(initCoord(t).allCoord(goodSpotIdx,1:3));
            % count fit
            ct = ct + 1;
            % remember current bad spots
            badSpotIdxLOld(:,ct) = badSpotIdxL;

        end

    end
end % loop good frames

% loop to get the "between frames" - stuff
for t=[metaphaseFrames(1:end-1),metaphaseFrames(2:end)]'
    % p-value of distribution comparison
    [dummy,planeFit(t(1)).deltaP] = kstest2(...
        planeFit(t(1)).planeCoord(planeFit(t(1)).inlierIdx,1),...
        planeFit(t(2)).planeCoord(planeFit(t(2)).inlierIdx,1));

    % change in plane orientation
    planeFit(t(1)).deltaAngle = acos(dot(planeFit(t(1)).planeVectors(:,1),...
        planeFit(t(2)).planeVectors(:,1))) *180/pi;
end

% plot everything

% eigenRatio
if verbose > 0
    figure('Name',sprintf('EigenRatio and distribution for %s',dataStruct.projectName))
    time = (1:nTimepoints)'.*dataStruct.dataProperties.timeLapse/60; % time in minutes
    subplot(2,1,1)
    plot(time,10*eigenValues(:,1)./mean(eigenValues(:,2:3),2),'-ro',...
        time,eigenValues(:,1),'--m',...
        time,mean(eigenValues(:,2:3),2),'--b+',...
        time,eigenValues(:,2),'.b',...
        time,eigenValues(:,3),'.b');
    legend('10x eigenRatio','different eigenValue','meanSimilar','similar1','similar2')
    xMinMax = [0,(nTimepoints+1)*dataStruct.dataProperties.timeLapse/60];
    xlim(xMinMax)
    ylim([0 50])
    xlabel('Time (min)')
    hold on
    % plot lines showing cutoff and mitosis
    plot(xMinMax,[10 10],':k')
    plot(metaphaseFrames([1 1])*dataStruct.dataProperties.timeLapse/60,[0 20],'k',...
        metaphaseFrames([end end])*dataStruct.dataProperties.timeLapse/60,[0 20],'k')
    text(metaphaseFrames(1)*dataStruct.dataProperties.timeLapse/60+1,13,'Mitosis')

    % distribution parameters
    distParms = cat(3,planeFit.distParms);
    subplot(2,1,2)
    plot(repmat(time(metaphaseFrames),1,3),squeeze(distParms(1:3,1,:))')
    hold on
    plot(time(metaphaseFrames),10*squeeze(distParms(4,1,:))','k')
    xlim(xMinMax)
    legend('variance','skewness','kurtosis','10*p(normal)')

    figure('Name',sprintf('P-value and orientation changes for %s',dataStruct.projectName))
    time = (1.5:nTimepoints-0.5)'.*dataStruct.dataProperties.timeLapse/60; % time in minutes
    subplot(2,1,1)
    plot(time(metaphaseFrames(1:end-1)),-log10(cat(1,planeFit.deltaP)))
    title('p-value for change in distribution (-log(p))')
    % plot 5% and 0.1% thresholds
    hold on
    plot(xMinMax,-log10([0.05 0.05]),':k')
    text(time(1),-log10(0.05)+0.15,'5% threshold')
    plot(xMinMax,-log10([0.001 0.001]),':k')
    text(time(1),-log10(0.001)+0.15,'0.1% threshold')
    xlim(xMinMax)
    ylim([0,-log10(0.001)+0.5])
    subplot(2,1,2)
    plot(time(metaphaseFrames(1:end-1)),(cat(1,planeFit.deltaAngle)))
    title('orientation change (deg)')
    xlim(xMinMax)

    % histogram of spot distribution
    c=-5.75:0.5:5.75;
    nBins = length(c);
    z = zeros(nBins,nTimepoints);
    for t=metaphaseFrames',z(:,t)=hist(planeFit(t).planeCoord(:,1),c);end
    uiViewPanel,imshow(z,[]),cm=isomorphicColormap('b');colormap(cm)
end


% assign out
dataStruct.planeFit = planeFit;

warning(warningState);









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SUBFUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [eigenVectors, eigenValues, meanCoord] = eigenCalc(coordinates)
% find eigenVectors, eigenValues, meanCoordinates for plane fitting
coordCov = cov(coordinates);

meanCoord = mean(coordinates,1);

% eigenVectors/eigenValues
[eigenVectors,eVals] = eig(coordCov);
eigenValues = diag(eVals)';

% compare eigenValues
diffs = pdist(eigenValues');
% indices to understand lowest
[u,v] = find(tril(ones(3),-1));
[dummy,idx] = min(diffs);
% find two close, one far
closestIdx = [u(idx),v(idx)];
farIdx = setdiff(1:3,closestIdx);


% sort eigenVectors, eigenValues. X-component of first eigenvector should
% always be positive
eigenVectors = eigenVectors(:,[farIdx,closestIdx]).*sign(eigenVectors(1,farIdx));
eigenValues = eigenValues([farIdx,closestIdx]);



