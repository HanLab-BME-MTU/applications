function [haloCoord, extraNucCoord, imgHalo]= ptFindHalos (inputImage, haloMinSize, haloMaxSize)
% ptFindHalos detects light areas and tries to fit calls into them
%
% SYNOPSIS       [haloCoord,imgHalo] = ptFindHalos (inputImage, haloMinSize, haloMaxSize)
%
% INPUT          inputImage    : either original image or segmented image 
%                haloMinSize : minimal size for halo
%                haloMaxSize : maximal size for halo
%
% OUTPUT         haloCoord : coordinates that have been found
%                extraNucCoord: coordinates of halos that are really nuclei as well
%                imgHalo   : binary image giving the areas of our halos
%
% DEPENDENCIES   ptFindHalos uses { nothing }
%                                  
%                ptFindHalos is used by { ptTrackCells
%                                         ptInitializeJob }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Colin Glass           Feb 04          Initial release
% Andre Kerstens        Mar 04          Cleaned up source

% The input image is already segmented, so take the pixels with value 3 
imgHalo = inputImage == 3;

% Let's do some morphological ops now to clean the image up
imgHalo = imerode (imgHalo, strel ('disk', 2));
imgHalo = imdilate (imgHalo, strel ('disk', 2));
imgHalo = imopen (imgHalo, strel ('disk', 3));
imgHalo = imclose (imgHalo, strel ('disk', 2));

% Label groups
imgHaloLabeled = bwlabel (imgHalo);

% Calculate halo properties using the regionprops function
haloProperties = regionprops (imgHaloLabeled, 'Area', 'PixelList');

% Find really small and big regions and temporarily store this info
haloArea = [haloProperties.Area];
haloMinMaxArea = find ((haloArea < haloMinSize) | (haloArea > haloMaxSize*2));
tempProperties = haloProperties (haloMinMaxArea);

% Get all pixels within those regions
minMaxCoord = cat (1, tempProperties(:).PixelList);

% Set all of those pixels to zero so that we can forget about them
for iCount = 1 : size (minMaxCoord, 1)
   imgHalo (minMaxCoord(iCount, 2), minMaxCoord(iCount, 1)) = 0;
end

% Now label the image again so we can get all the coordinates out
imgHaloLabeled = bwlabel (imgHalo);

% And get the properties
haloProperties = regionprops (imgHaloLabeled, 'Centroid', 'Eccentricity', 'MajorAxisLength');

% Finally we store the coordinates of all the halos
haloCoord = round (cat (1, haloProperties.Centroid));

% Also get the eccentricity of the regions so that we can distinguish big
% round areas that are probably cell nuclei as well (eccentricity of 0 is
% a perfect circle). eccentricity of exactly 0 should be discarded.
haloEccentricity = [haloProperties.Eccentricity];
haloMajorAxisLength = [haloProperties.MajorAxisLength];
haloDiameter = 2.0 * sqrt (haloMaxSize / pi);
haloSmallEcc = find ((haloEccentricity < 0.5) & (haloEccentricity ~= 0) & ...
                     (haloMajorAxisLength < haloDiameter));
tempProperties = haloProperties (haloSmallEcc);
extraNucCoord = round (cat (1, tempProperties(:).Centroid));

% Really big areas also are likely to be nuclei, so let's find these now
% First some morphological operations 
%diskSize = round (sqrt (haloMaxSize / 10) / 2);
%imgHaloTemp = imclose (imgHalo, strel ('disk', 3));
%imgHaloTemp = imerode (imgHaloTemp, strel ('disk', diskSize));
%imgHaloTemp = imopen (imgHaloTemp, strel ('disk', diskSize));

% Label the resulting image again
%imgHaloTempLabeled = bwlabel (imgHaloTemp);

% Get the area properties
%haloProps = regionprops (imgHaloTempLabeled, 'Area', 'Centroid');

% Get rid of the really small areas in case these still exist and only keep
% the bigger ones
%if ~isempty (haloProps)
%   haloAreaTemp = [haloProps.Area];
%   haloMinAreaTemp = find (haloAreaTemp > haloMinSize);
%   tempProps = haloProps (haloMinAreaTemp);
%
%   % At last find the coordinates of these big areas
%   bigHaloCoord = round (cat (1, extraNucCoord, tempProps(:).Centroid));
%   %extraNucCoord = round (cat (1, extraNucCoord, bigHaloCoord));
%   extraNucCoord = bigHaloCoord;
%end
