function [nucCoord, imgNuclei] = ptFindNuclei (inputImage, nucleiMinSize, nucleiMaxSize)
% ptFindNuclei detects dark areas and tries to fit cells into them
%
% SYNOPSIS       [nucCoord, imgNuclei] = ptFindNuclei (inputImage, nucloiMinSize, nucloiMaxSize)
%
% INPUT          inputImage    : either original image or segmented image (depends on method)
%                nucleiMinSize : minimal size for nuclei
%                nucleiMaxSize : maximal size for nuclei
%
% OUTPUT         nucCoord         : coordinates found
%                imgNuclei : binary image showing the areas of nuclei
%
% DEPENDENCIES   ptFindNuclei uses { nothing }
%                                  
%                ptFindNuclei is used by { ptTrackCells
%                                               ptInitializeJob } 
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Colin Glass           Feb 04          Initial release
% Andre Kerstens        Mar 04          Cleaned up source

imgNuclei = inputImage == 1;

% Let's do some morphological ops now to clean the image up
imgNuclei = imerode (imgNuclei, strel ('disk', 2));
imgNuclei = imdilate (imgNuclei, strel ('disk', 2));
imgNuclei = imopen (imgNuclei, strel ('disk', 3));
imgNuclei = imclose (imgNuclei, strel ('disk', 3));

% Label the image
imgNucleiLabeled = bwlabel (imgNuclei);

% Calculate nucloi properties using the regionprops function from the image toolbox
nucleiProperties = regionprops (imgNucleiLabeled, 'Area', 'PixelList', 'Eccentricity');

% Find really small and big regions and temporarily store this info
nucleiArea = [nucleiProperties.Area];
nucleiMinMaxArea = find ((nucleiArea < nucleiMinSize) | (nucleiArea > nucleiMaxSize*2));
tempProperties = nucleiProperties (nucleiMinMaxArea);

% Get all pixels within those regions
minMaxCoord = cat (1, tempProperties(:).PixelList);

% Also get the eccentricity of the regions so that we can remove long areas that are not
% nuclei, but stuff that comes out of the nuclei. 0.98 should be a good value to find these.
nucleiEccentricity = [nucleiProperties.Eccentricity];
nucleiBigEcc = find (nucleiEccentricity > 0.98);
tempProperties2 = nucleiProperties (nucleiBigEcc);

% Get the pixels in those areas as well
bigEccCoord = cat (1, tempProperties2(:).PixelList);

% Concatenate all of the coordinates to be removed together
removeCoord = cat (1, minMaxCoord, bigEccCoord);

% Set all of those pixels to zero so that we can forget about them
for iCount = 1 : size (removeCoord, 1)
   imgNuclei (removeCoord(iCount, 2), removeCoord(iCount, 1)) = 0;
end

% Now label the new image again
imgNucleiLabeled = bwlabel (imgNuclei);

% Calculate nuclei properties using the regionprops function again
nucleiProperties = regionprops (imgNucleiLabeled, 'Centroid');

% And get the coordinates of all the remaining nuclei
nucCoord = round (cat (1, nucleiProperties.Centroid));
