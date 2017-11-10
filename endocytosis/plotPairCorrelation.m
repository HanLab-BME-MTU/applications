function [experiment] = plotPairCorrelation(experiment,dist,plotMask)

% plotPairCorrelation calculates the density of pits defined by rest
% and that fall within an inputMask
%
% INPUT:   experiment=   data structure pointing to all the
%                       image/detection/tracking data for a given condition;
%                       for e.g. 10 cell movies, the structure should have 10
%                       entries; the data structure can be created with the
%                       function loadConditionData, and it needs to contain
%                       at least the field
%                       .source, which is the (path) location of the
%                       lifetime information folder
%                       .framerate, which is the movie framerate, which is
%                       necessary for the lifetime restriction
%                       .status (optional), which contains a number for
%                       each pit which identifies it as belonging to a
%                       certain population, so that all pits that belong to
%                       population1 have a status value 1 (this is useful
%                       in calculating pair correlations for the different
%                       populations)
%                       [... minint maxint minmot maxmot]
%           dist =      [optional] distance vector for binning (default is
%                       1:20)
%           plotMask =  [optional] 1 to plot area mask and detections 0 to
%                       not (default is 0)
%
% OUTPUT
%           .pairCorrelation =  pair correlation as column vecotrs for each
%           movie; if pairCorrelation is a matrix then each column
%           corresponds to a given population defined by status vector and
%           sorted by value of status so that the first column is the pair
%           correlation for the lowest group number
% Uses:
%       makeCellMaskDetections
%       RipleysKfunction
%       calculatePitDenFromLR
%       makeCorrFactorMatrix
%
% Daniel Nunez, updated May 01, 2013

%%
%AREA MASK PARAMETERS
closureRadius = 20;
dilationRadius = 5;
doFill = 0;

%%
%interpret inputs
if nargin < 2 || isempty(dist)
    dist = 1:20;
end
if nargin < 3 || isempty(plotMask)
    plotMask = 0;
end


for iexp = 1:length(experiment)
    
    % framerate
    framerate = experiment(iexp).framerate;
    % image size
    imsize  = experiment(iexp).imagesize;
    
    tracks = loadTracks(experiment(iexp));
    % all positions
    maskPositions = [];
            maskPositions(:,1) = cell2mat(arrayfun(@(var)var.x(1,:),tracks,'UniformOutput',false))';
            maskPositions(:,2) = cell2mat(arrayfun(@(var)var.y(1,:),tracks,'UniformOutput',false))';
    %all track nucleations
    mpm = nan(length(tracks),2);
            mpm(:,1) = cell2mat(arrayfun(@(var)var.x(1,1),tracks,'UniformOutput',false))';
            mpm(:,2) = cell2mat(arrayfun(@(var)var.y(1,1),tracks,'UniformOutput',false))';
    
    %MAKE MASK
    imsizS = [imsize(2) imsize(1)];
    [areamask] = makeCellMaskDetections(maskPositions,closureRadius,...
        dilationRadius,doFill,imsize,plotMask,[]);
    %CALCULATE NORMALIZED AREA FROM MASK
    normArea = bwarea(areamask);
    % CREATE CORRECTION FACTOR MATRIX FOR THIS MOVIE using all objects
    corrFacMat = makeCorrFactorMatrix(imsizS, dist, 10, areamask');
    
    %CALCULATE PIT DENSITY
    %if there is a status field used to divide pit populations
    %     if isfield(experiment(iexp),'status')
    %         %for each value in the status vector
    %         pair = nan(length(dist),length(unique(experiment(iexp).status)));
    %         populations = unique(experiment(iexp).status);
    %         for ipop = 1:length(populations)
    %             mpm1 = [full(matX(findPos & repmat(experiment(iexp).status == populations(ipop),size(matX,2),1)'))...
    %                 full(matY(findPos & repmat(experiment(iexp).status == populations(ipop),size(matX,2),1)'))];
    %             [kr,lr]=RipleysKfunction(mpm1,mpm1,imsizS,dist,corrFacMat,normArea);
    %             [currDen] = calculatePitDenFromLR(kr,dist);
    %             pair(:,ipop) = currDen;
    %         end
    %     else
    %
    
    [kr,lr,pair]=RipleysKfunction(mpm,mpm,imsizS,dist,corrFacMat,normArea);
    
    %store pair in each experiment structure
    experiment(iexp).pairCorrelation = pair;
    experiment(iexp).clustering = sum(pair(1:2));
    experiment(iexp).nucleationDensity = length(tracks)/normArea; 
    
    %SCRAMBLE MPM
    [mpm] = makeRandomMPM(mpm, areamask',1);
    [kr,lr,pair]=RipleysKfunction(mpm,mpm,imsizS,dist,corrFacMat,normArea);
    experiment(iexp).pairCorrelationRandom = pair;
    
    
end

end %of function