function [windowInfo] = plotPairCorrelation_SpatialWindows(experiment,dist,rest,plotMask,windowPolygons);

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
%           rest    =   [optional] restriction vector can have variable length;
%                       minimum length is five, where the entries are
%                       [stat da minfr minlft maxlft]
%                       (default is [1 1 4 1 300])
%                       optionally, the length can be extended to nine,
%                       where the additional entries are
%                       [... minint maxint minmot maxmot]
%           dist =      [optional] distance vector for binning (default is
%                       1:20)
%           plotMask =  [optional] 1 to plot area mask and detections 0 to
%                       not (default is 0)
%
% OUTPUT
%           pairCorr =  pair correlation as column vecotrs for each movie
%           clustering = sum of pairCorrelation for the first two pixels
%                       for ech movie
% Uses:
%       makeCellMaskDetections
%       RipleysKfunction
%       calculatePitDenFromLR
%       makeCorrFactorMatrix
%
% Daniel Nunez, updated April 24, 2013

%% EXPLANATION of restriction values:
% rest = [stat da minfr minlft maxlft minint maxint minmot maxmot]
%
% stat  =   object status, possible values are 1, 2, 3
%           1: object appears and disappears in the movie
%           2: object is present for entire length of the movie
%           3: object is present in either first or last frame of the movie
% da    =   appearance/disappearance status, possible values are 1,0,-1
%           1: marks first frame
%           0: marks middle frame
%           -1: marks last frame
% minfr =   minimum lifetime/length in FRAMES - e.g. 4, to exclude tracking
%           artifacts of false detection positives
% minlft =  minimum lifetime in SECONDS - e.g. 60 to select for productive
%           clathrin-coated pits
% maxlft =  maximum lifetime in SECONDS - e.g. 25 to select for abortive
%           clathrin-coated pits
%
% OPTIONAL:
%
% minint =  minimum normalized intensity (ranging from 0 to 1)
% maxint =  maximum normalized intensity (ranging from 0 to 1)
% minmot =  minimum normalized motility (ranging from 0 to 1)
% maxmot =  maximum normalized motility (ranging from 0 to 1)
%
% These latter criteria allow you to select e.g. the brightest 10% of the
% population, or the faster 50%.
%
% EXAMPLE:  to select the positions where productive pits appear
%           rest1 = [1 1 4 60 300]
%           to select the positions where abortive pits are located in each
%           frame
%           rest2 = [1 0 4 8 25]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%AREA MASK PARAMETERS
closureRadius = 20;
dilationRadius = 5;
doFill = 0;

%%
%interpret inputs
if nargin < 2 || isempty(dist)
    dist = 1:20;
end
if nargin < 3 || isempty(rest)
    rest = [1 1 4 1 300];
end
if nargin < 4 || isempty(plotMask)
    plotMask = 0;
end
if nargin < 5 || isempty(windowPolygons)
    windowPolygons = [];
end


for iexp = 1:length(experiment)
    
    tracks = loadTracks(experiment(iexp));
    imsizS = [experiment(iexp).imagesize(2) experiment(iexp).imagesize(1)];
    
    %count windows
    windowsPerSlice = cellfun(@(x)numel(x),windowPolygons);
    
    %pick out pits that fall in windows
    for i = 1:length(windowsPerSlice)
        for j = 1:windowsPerSlice(i)
            
            mpm = nan(length(tracks),2);
            mpm(:,1) = cell2mat(arrayfun(@(var)var.x(1,1),tracks,'UniformOutput',false))';
            mpm(:,2) = cell2mat(arrayfun(@(var)var.y(1,1),tracks,'UniformOutput',false))';
            
            border = [windowPolygons{i}{j}{:}];
            in  = inpolygon(mpm(:,1),mpm(:,2),border(1,:),border(2,:));
            
            windowInfo(i,j).pits = mpm(in,:);
            
            if length(windowInfo(i,j).pits) >= 3
                %MAKE MASK
                [areamask] = makeCellMaskDetections(windowInfo(i,j).pits,closureRadius,dilationRadius,doFill,experiment(iexp).imagesize,0,[]);
                %CALCULATE NORMALIZED AREA FROM MASK
                normArea = bwarea(areamask);
                
                % CREATE CORRECTION FACTOR MATRIX FOR THIS MOVIE using all objects
                corrFacMat = makeCorrFactorMatrix(imsizS, dist, 10, areamask');
                
                %CALCULATE PIT DENSITY
                [kr,lr,currDen]=RipleysKfunction(windowInfo(i,j).pits,windowInfo(i,j).pits,imsizS,dist,corrFacMat,normArea);
                
                
                %calculate other good stuff
                windowInfo(i,j).pairCorrelation = currDen;
                windowInfo(i,j).nucDen = size(windowInfo(i,j).pits,1)/normArea;
                windowInfo(i,j).lft = [tracks(in).lifetime_s];
            else
                windowInfo(i,j).pairCorrelation = [];
                windowInfo(i,j).nucDen = [];
                windowInfo(i,j).lft = [];
            end
        end
        
    end
end %of function


