function [experiment] = determineInitiationDensity(experiment,rest,plotMask,inputMask)

% determineInitiationDensity calculates the density of pits defined by rest
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
%           inputMask = binary mask; pits that fall within pixels of value
%                       one will be counted towards density.
% OUTPUT
%           experiment.initiationDen in pits/frame/pixel
%           experiment.initiationDen in pits/second/micrometer^2
%
% Uses:
%       makeCellMaskDetections
%
% Daniel Nunez, updated May 5, 2009

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

%%
%AREA MASK PARAMETERS
closureRadius = 20;
dilationRadius = 5;
doFill = 0;

%%
%interpret inputs
if nargin < 2 || isempty(rest)
    rest = [1 1 4 1 300];
end
if nargin < 3 || isempty(plotMask)
   plotMask = 0;
end
if nargin < 4 || isempty(inputMask)
    inputMask = [];
end

%convert movie paths to correct OS
    %[experiment]=changePathUsingEndocytosis(experiment);


%make image binary
inputMask = im2bw(inputMask);

for iexp = 1:length(experiment)

    %Load Lifetime Information
    load([experiment(iexp).source filesep 'LifetimeInfo' filesep 'lftInfo.mat']);
    % status matrix
    statMat = lftInfo.Mat_status;
    % lifetime matrix
    lftMat = lftInfo.Mat_lifetime;
    % x-coordinate matrix
    matX = lftInfo.Mat_xcoord;
    % y-coordinate matrix
    matY = lftInfo.Mat_ycoord;
    % disapp status matrix
    daMat = lftInfo.Mat_disapp;
    % framerate
    framerate = experiment(iexp).framerate;
    %movie length
    movieLength = experiment(iexp).movieLength;
    imSize = experiment(iexp).imagesize;
    %use status to select for a population of pits
    if isfield(experiment,'status')
        status = experiment(iexp).status;
    else
        status = ones(1,size(daMat,1));
    end
    
    %find all pits in movie that meet requirements specified by restriction
    %vector
    findPos = find((statMat==rest(1,1)) & (daMat==rest(1,2)) &...
        (lftMat>rest(1,3)) & (lftMat>round(rest(1,4)/framerate)) & (lftMat<round(rest(1,5)/framerate)) &...
        repmat(status',1,size(statMat,2)) == 1);
    

    %get pits inside mask
    if exist('inputMask','var') && ~isempty(inputMask)
        findPos = findPos(diag(inputMask(ceil(matY(findPos)),ceil(matX(findPos)))) == 1);
    end

    [areamask] = makeCellMaskDetections([matX(~isnan(matX)),matY(~isnan(matY))],closureRadius,dilationRadius,doFill,imSize,plotMask,[]);
    normArea = bwarea(areamask);
    
    experiment(iexp).initiationDen = length(findPos)/normArea/movieLength;
    experiment(iexp).initiationDenUnits = length(findPos)/(normArea*0.067^2)/(movieLength*framerate)*60;

end %for each experiment

end %of function