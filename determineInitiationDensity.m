function [experiment] = determineInitiationDensity(experiment,rest,inputMask);

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
%           rest    =   restriction vector can have variable length;
%                       minimum length is five, where the entries are
%                       [stat da minfr minlft maxlft]
%                       optionally, the length can be extended to nine,
%                       where the additional entries are
%                       [... minint maxint minmot maxmot]
%           inputMask = binary mask; pits that fall within pixels of value
%           one will be counted towards density.
% OUTPUT
%           experiment.initiationDen in pits/frame/pixel
%           experiment.initiationDen in pits/second/micrometer^2
%
% Uses:
%       determineMovieLength
%       determineImagesize
%       makeCellMaskDetections
%
% Daniel Nunez, updated March 11, 2009

%interpret results
%make image binary
inputMask = im2bw(inputMask);

%Fill in Missing Data
%needed to normalize density by movie length
[experiment] = determineMovieLength(experiment);
%needed to make mask
[experiment] = determineImagesize(experiment);

for iexp = 1:length(experiment)
    
    %gives wait bar; comment if not wanted
    waitHandle = waitbar(iexp/length(experiment),['running movie ' num2str(iexp) ' out of ' num2str(length(experiment))]);
    
    %Load Lifetime Information
    cd([experiment(iexp).source filesep 'LifetimeInfo'])
    load('lftInfo')
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
    
    %find all pits in movie that meet requirements specified by restriction
    %vector
    findPos = find((statMat==rest(1,1))& (daMat==rest(1,2)) &...
        (lftMat>rest(1,3)) & (lftMat>round(rest(1,4)/framerate)) & (lftMat<round(rest(1,5)/framerate)));

    %get pits inside mask
    if exist('inputMask','var') && ~isempty(inputMask)
    findPos = findPos(diag(inputMask(ceil(matY(findPos)),ceil(matX(findPos)))) == 1);
    end
    
    [areamask] = makeCellMaskDetections([matX(findPos),matY(findPos)],40,imSize);
    normArea = bwarea(areamask);  
    
    experiment(iexp).initiationDen = length(findPos)/normArea/movieLength;
    experiment(iexp).initiationDenUnits = length(findPos)/(normArea*0.067^2)/(movieLength*framerate)*60;
    
    close(waitHandle)
    
end %for each experiment

end %of function