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

% Close gaps which are smaller than a disk with radius 3
%imgNuclei = imclose(imgNuclei, strel('disk',3));

% And dilate the picture as well (same disk size)
%imgNuclei = imdilate(imgNuclei, strel('disk',3));

% Label results of the processed image
%imgNucleiLabeledBefore = bwlabel (imgNuclei);
%lastNucloiGroup = max (max (imgNucleiLabeledBefore)) ;

% More morphological operations, but this time the regions are already labeled.
% So if one region gets seperated thus becoming two or more, it still has the same numbers.
%imgNuclei = imerode (imgNuclei, strel('disk',7));
%imgNuclei = imdilate (imgNuclei, strel('disk',6));

% Now we label the image again. So we have one image labeled before the
% second morph and one labeled after
%imgNucleiLabeledAfter = bwlabel (imgNuclei);

%counter = 0;
%membership = [];

% For every group of imgNucleiLabeledBefore (labeled before second morph oper) we allow
% one labeled group in imgNucleiLabeledAfter (labeled after...). If there are more then 
% one region we selected the biggest.
% for group = 1 : lastNucloiGroup
%    counter = counter + 1;
%    nucloiGroup = imgNucleiLabeledAfter (find (imgNucleiLabeledBefore == group));
%    nucloiGroup = nucloiGroup (find (nucloiGroup));
%    uniqueNucloiGroup = unique (nucloiGroup);
%    if ~isempty (uniqueNucloiGroup)
%       if length (uniqueNucloiGroup) > 1.1  
%          [uniqueEntries, numberOfOccurences] = countEntries (nucloiGroup);
%          [dummy, groupIndex] = max (numberOfOccurences);
%          membership (group, 1) = uniqueEntries (groupIndex);
%       else 
%          membership (counter, 1) = uniqueNucloiGroup;
%       end
%    else
%       counter = counter - 1;
%    end
% end

% In membership we have at most one group (belonging to imgNucleiLabeledAfter) for
% every each group belonging to imgNucleiLabeledBefore
% imgNuclei is overwritten now only including the groups (of imgNucleiLabeledAfter) given in membership
% if ~isempty (membership)
%    % find all the non-zero groups
%    membership = membership (find (membership));
%    imgNuclei = ismember (imgNucleiLabeledAfter, membership);
% 
%    % Label the new image again
%    clear imgNucleiLabeled;
%    imgNucleiLabeled = bwlabel (imgNuclei);
% end

