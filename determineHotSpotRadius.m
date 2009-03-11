function [hotSpotRadius] = determineHotSpotRadius(experiment,rest);

% determineHotSpotRadius calculates radius of hot spot based on pair
% correlation (takes mean value of pair correlation from pixels 10 to 20
% and finds the radius at which the pair correlation is equal to this value
% plus one)
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
%                       
% OUTPUT
%           hotSpotRadius: mean value of pair correlation from pixels 10 to   
%           20 plus one)
%
% Uses: determineMovieLength
%       determineImagesize
%       RipleysKfunction
%
% Daniel Nunez, January 20, 2009

dist = 1:1:20;

%Fill in Missing Data
[experiment] = determineMovieLength(experiment);
[experiment] = determineImagesize(experiment);

for iexp = 1:length(experiment)
    
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
    % image size
    imsize  = experiment(iexp).imagesize;
    
    %find all pits in movie that meet requirements specified by restriction
    %vector
    findPos = find((statMat==rest(1,1)) & (daMat==rest(1,2)) &...
        (lftMat>rest(1,3)) & (lftMat>round(rest(1,4)/framerate)) & (lftMat<round(rest(1,5)/framerate)));
    
    msx = imsize(1);
    msy = imsize(2);
    imsizS = [imsize(2) imsize(1)];
    % construct convex hull out of complete point distribution
    % combined mpm
    selx = full(matX); selx = selx(isfinite(selx)); selx = nonzeros(selx(:));
    sely = full(matY); sely = sely(isfinite(sely)); sely = nonzeros(sely(:));
    combMPM = [ selx sely ];
    K = convhull(combMPM(:,1),combMPM(:,2));
    % edge points of the convex hull
    cpointsx = combMPM(K,1);
    cpointsy = combMPM(K,2);
    % create mask
    areamask = poly2mask(cpointsx,cpointsy,msx,msy);
    % CREATE CORRECTION FACTOR MATRIX FOR THIS MOVIE using all objects
    normArea = sum(areamask(:));

    mpmPits = [full(matX(findPos)) full(matY(findPos))];
    [kr,lr]=RipleysKfunction(mpmPits,mpmPits,imsizS,dist,[],normArea);
    Lpits(:,iexp) = lr;  
end

    [dlx,dly,dlz] = size(Lpits);
    carea = dist.^2;
    careadiff = carea; careadiff(2:length(careadiff)) = diff(carea);
    amat = repmat(careadiff',1,dly);
    dmat = repmat(dist',1,dly);
    currLR = Lpits;
    currKR = (currLR+dmat).^2;
    currKRdiff = currKR;
    currKRdiff(2:length(dist),:) = diff(currKR,1);
    currDen = currKRdiff./amat;
    pitDen = nanmedian(currDen,2);

%fit spline to productive data
pp = csaps(dist,pitDen);
denFitValues = fnplt(pp);
%find max value of fit
maxNormDen = max(denFitValues(2,:));
%find density value where radius is 10 pixels
findMidRad = find(denFitValues(1,:) >= 10,1,'first');
%take the mean density from 10 pixels to the end as the far distance value
farDen = mean(denFitValues(2,findMidRad:end));
%get point at 10 percent of difference in between maximum and long distance
%densities
findRad = find(denFitValues(2,:) < 1+farDen,1,'first');
hotSpotRadius = denFitValues(1,findRad);

end %of function