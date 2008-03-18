function makiEvaluateCondition(cutoff,serverType)
%MAKIEVALUATECONDITION evaluates the three-condition movies to find good imaging parameters
%
% SYNOPSIS: makiEvaluateCondition(cutoff)
%
% INPUT :   Cutoff(opt) : displacement cutoff for point correspondences
%           serverType  : ID to indicate gui-loading path
%
%           The function will load three movies with frames taken in
%           the order 123123123...
%
% OUTPUT
%
% REMARKS makiEvaluateCondition returns the percentage of false positives
%           and false negatives in the initial spot detection using the
%           movie with the highest exposure (movie number three) as ground
%           truth
%
% created with MATLAB ver.: 7.4.0.287 (R2007a) on Windows_NT
%
% created by: jdorn
% DATE: 03-Jul-2007
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def_cutoff = 1; % cutoff in microns
if nargin < 1 || isempty(cutoff)
    cutoff = def_cutoff;
end
if nargin < 2
    serverType = [];
end

% load data files
%done = false;
data = struct('dataStruct',[]);
for i = 1:3
    h = warndlg(sprintf('Please load movie %i. Longest exposure should be last!',i));
    uiwait(h);
    data(i).dataStruct = makiLoadDataFile(serverType);
end

% for every frame: count spots. Also, LAP the frames onto the movie with
% the longest exposure and count false positives (points not in ground
% truth) and false negatives (points only in ground truth). Maximum
% displacement allowed is 500 nm

pix2mu = [data(1).dataStruct.dataProperties.PIXELSIZE_XY,data(1).dataStruct.dataProperties.PIXELSIZE_XY,...
    data(1).dataStruct.dataProperties.PIXELSIZE_Z];
imSize = [data(1).dataStruct.movieHeader.numRows,data(1).dataStruct.movieHeader.numCols,data(1).dataStruct.movieHeader.numZSlices];

% all dataProperties are the same
nData = 3;%length(data);
for iData = nData:-1:1; % count backwards so that we have nSpots defined on time
    % get nSpots
    dataStruct = data(iData).dataStruct;
    nTimepoints = length(dataStruct.initCoord);
    data(iData).nSpots = [...
        catStruct(1,'dataStruct.initCoord.nSpots'),...
        zeros(nTimepoints,2)]; %number, false pos, false neg, 100%

    % check for cropping
    if isempty(dataStruct.dataProperties.crop)
        dataStruct.dataProperties.crop = zeros(2,3);
    end
    crop = dataStruct.dataProperties.crop(:,1:3);
    isCrop = any(crop,1);
    crop(1,~isCrop) = 1;
    crop(2,~isCrop) = imSize(~isCrop);

    for t = 1:nTimepoints



        % compare coordinates to #3
        if iData == 3
            coords = dataStruct.initCoord(t).allCoordPix(:,1:3);
            coordTruth = coords + repmat(crop(1,:)-1,size(coords,1),1);
            coordTruth = coordTruth.*repmat(pix2mu,size(coords,1),1);
        else


            coords = dataStruct.initCoord(t).allCoordPix(:,1:3);
            coordTest = coords + repmat(crop(1,:)-1,size(coords,1),1);
            coordTest = coordTest.*repmat(pix2mu,size(coords,1),1);



            % get distance
            dm = createDistanceMatrix(coordTest,coordTruth);
            % cutoff at 1 um
            dm(dm>cutoff) = -1;
            % lap
            [toTruth,fromTruth] = lap(dm, -1, [], 1);
            % false positives: not linked to truth
            data(iData).nSpots(t,2) = sum(toTruth(1:data(iData).nSpots(t,1))>data(3).nSpots(t,1));
            % false negatives: not linked from truth
            data(iData).nSpots(t,3) = sum(fromTruth(1:data(3).nSpots(t,1))>data(iData).nSpots(t,1));
            % 100%: nSpots of ground truth
            data(iData).nSpots(t,4) = data(3).nSpots(t,1);

        end
    end
end

% figure,
% plot(data(1).nSpots,'xr'),
% hold on,
% plot(data(2).nSpots,'+g'),
% plot(data(3).nSpots,'ob')

% evaluate: plot percent fp/fn
figure('Name',sprintf('%s cutoff %1.1f',data(3).dataStruct.projectName,cutoff));
for i=1:2
    subplot(2,2,i)
    plot(1:nTimepoints,100*(data(i).nSpots(:,2))./data(i).nSpots(:,4),'+r',...
        1:nTimepoints,100*(data(i).nSpots(:,3))./data(i).nSpots(:,4),'ob');
    ylim([0,20]);
    hold on
    avgNeg = 100*mean(data(i).nSpots(:,3)./data(i).nSpots(:,4));
    avgPos = 100*mean(data(i).nSpots(:,2)./data(i).nSpots(:,4));
    plot([1,nTimepoints],[avgPos,avgPos],'r',...
        [1,nTimepoints],[avgNeg,avgNeg],'b')
    title(data(i).dataStruct.projectName)
    legend('false positives%','false negatives%')
end
subplot(2,2,3),plot(data(1).nSpots(:,1),'xr'),
hold on,
plot(data(2).nSpots(:,1),'+b'),
plot(data(3).nSpots(:,1),'og')
title('nSpots')
legend(data(1).dataStruct.projectName,data(2).dataStruct.projectName,data(3).dataStruct.projectName)




