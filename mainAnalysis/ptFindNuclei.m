function [coord, imgNuclei] = ptFindNuclei (inputImage, nucloiLevel, nucloiMinSize, nucloiMaxSize, method)
% ptFindNuclei detects dark areas and tries to fit cells into them
%
% SYNOPSIS       [coord, imgNuclei] = ptFindNuclei (inputImage, nucloiLevel, nucloiMinSize, nucloiMaxSize, method)
%
% INPUT          inputImage    : either original image or segmented image (depends on method)
%                nucloiLevel   : level used for minima detection
%                nucloiMinSize : minimal size for nuclei
%                nucloiMaxSize : maximal size for nuclei
%                method        : 1 or 2. Says if clustering or image segmentation has been applied 
%                                (changes what ptFindNuclei actually does)
%
% OUTPUT         coord         : coordinates found
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

% Do some input checking first to determine the segmentation method used
if method == 1          % Clustering
   % The input image is already segmented, so take the pixels with value 3
   imgNuclei = inputImage == 1;
    
elseif method == 2 || method == 3     % Thresholding
   % Calculate the binary image based on the nucloi intensity level
   imgNuclei = imextendedmin (inputImage, nucloiLevel);
    
else
   fprintf (1, 'ptFindNuclei: invalid method. Please use method 1 or 2.\n');
end 

% Let's do some morphological ops now to clean the image up
imgNuclei = imfill (imgNuclei,'holes');
imgNuclei = imerode (imgNuclei, strel ('disk', 5));
imgNuclei = imdilate (imgNuclei, strel ('disk', 5));
imgNuclei = bwareaopen (imgNuclei, 50);

% Label the image
imgNucleiLabeled = bwlabel (imgNuclei);

% Calculate nucloi properties using the regionprops function from the image toolbox
nucloiProperties = regionprops (imgNucleiLabeled, 'Area', 'PixelList');

% Find really small and big regions and temporarily store this info
nucleiArea = [nucloiProperties.Area];
nucloiMinMaxArea = find ((nucleiArea < round (nucloiMinSize / 4 * 3)) | (nucleiArea > nucloiMaxSize));
tempProperties = nucloiProperties (nucloiMinMaxArea);

% Get all pixels within those regions
minMaxCoord = cat (1, tempProperties(:).PixelList);

% Set all of those pixels to zero so that we can forget about them
for iCount = 1 : size (minMaxCoord, 1)
   imgNuclei (minMaxCoord(iCount, 2), minMaxCoord(iCount, 1)) = 0;
end

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

% Calculate nucloi properties using the regionprops function from the image toolbox
nucloiProperties = regionprops(imgNucleiLabeled, 'Area', 'PixelList', 'MajorAxisLength', 'Centroid', 'Eccentricity');

% Remove really small nuclei
clear nucleiArea; clear tempProperties;
nucleiArea = [nucloiProperties.Area];
nucleiMinArea = find (nucleiArea > nucloiMinSize);
tempProperties = nucloiProperties (nucleiMinArea);

% Remove really large nuclei
clear nucleiArea;
nucleiArea = [tempProperties.Area];
nucleiMaxArea = find(nucleiArea < nucloiMaxSize);
tempProperties2 = tempProperties (nucleiMaxArea);

coord = round (cat (1, tempProperties2.Centroid));
