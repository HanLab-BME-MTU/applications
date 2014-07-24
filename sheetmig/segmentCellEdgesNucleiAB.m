function [maskNuclei,watershedEdges,centersNuclei,centersCells] = ...
    segmentCellEdgesNuclei(imageN,imageC,plotRes)
%SEGMENTCELLEDGESNUCLEI segments nuclei by thresholding and cell edges by watershed
%
%SYNOPSIS [maskNuclei,watershedEdges,centersNuclei,centersCells] = ...
%    segmentCellEdgesNuclei(imageN,imageC,plotRes)
%
%INPUT  imageN    : Image of nuclei.
%       imageC    : Image of cell edges.
%       plotRes   : 1 to plot results, 0 otherwise.
%                   In the plot, nuclei are shown in blue, cell edges in
%                   green, the watershed edges in red, the cell centers as
%                   red pluses, the nuclear positions as yellow crosses.
%                   Optional. Default: 0.

%
%OUTPUT maskNuclei    : Mask of nuclei. 1 inside nuclei, 0 outside.
%       watershedEdges: Watershed of cell edge image. 0 indicates the
%                       edges, while the cell regions are numbered 1 ... N.
%       centersNuclei : X- and Y-coordinates of the nuclei.
%       centersCells  : X- and Y-coordinates of the cells.
%
%Khuloud Jaqaman, June 2010

% Input

%check whether images have been input
if nargin < 2
    error('--segmentCellEdgesNuclei: Please enter images of nuclei and cell edges.')
end

if nargin < 3 || isempty(plotRes)
    plotRes = 0;
end

% Nuclear segmentation by thresholding

%convert images to double
imageN = double(imageN);
imageC = double(imageC);

% filter nuclei image with a large Gaussian
imageNfilter = filterGauss2D(imageN, 10);


%get nuclear masks by thresholding
maskNuclei = double(blobSegmentThreshold(imageNfilter,500,plotRes));
maskNucleiLabel = bwlabel(maskNuclei);

%get the nuclear centers
STATS = regionprops(maskNucleiLabel,'Centroid'); %#ok<MRPBW>
centersNuclei = round(vertcat(STATS.Centroid));


% Cell edge segmentation by watershed

%filter cadherin image
imageCfilter = steerableFiltering(imageC,2,10);

%make image of nuclear centers to initialize the watershed
M = zeros(size(maskNuclei));
M(sub2ind(size(maskNuclei), centersNuclei(:,2), centersNuclei(:,1))) = 1;


%do the watershed
watershedEdges = watershedXT(imageCfilter-maskNuclei,4,M);

%get the cell centers
STATS = regionprops(watershedEdges,'Centroid'); 
centersCells = round(vertcat(STATS.Centroid));

% Display results

if plotRes
    
    Ldisp = watershedEdges;
    Ldisp(Ldisp>0) = -1;
    image3Color(:,:,3) = (imageN-min(imageN(:)))/(max(imageN(:))-min(imageN(:)));
    image3Color(:,:,2) = (imageC-min(imageC(:)))/(max(imageC(:))-min(imageC(:)));
    image3Color(:,:,1) = (Ldisp-min(Ldisp(:)))/(max(Ldisp(:))-min(Ldisp(:)));
    figure, imshow(image3Color,[])
    hold on
    plot(centersNuclei(:,1),centersNuclei(:,2),'yx')
    plot(centersCells(:,1),centersCells(:,2),'r+','MarkerSize',5)
    
end


% ~~~ the end ~~~


% Snippets of code from initial attempts

% %get region statistics
% STATS = regionprops(L,maskNucleiLabel,'MaxIntensity');
% nuclearInt = vertcat(STATS.MaxIntensity);
% 
% %merge objects that share a nuclear label
% intUnique = unique(nuclearInt);
% numUnique = length(intUnique);
% for i = 2 : numUnique
%     
%     indxCell = find(nuclearInt==intUnique(i));
%     numCells = length(indxCell);
%     
%     if numCells > 1
%        
%         tmpL = ismember(L,indxCell);
%         SE = strel('square',3);
%         tmpDilate = imdilate(tmpL,SE);
%         tmpErode = imerode(tmpDilate,SE);
%         
%         L(tmpErode) = indxCell(1);
%         
%     end
%     
% end

% %get threshold based on product of intensity and area
% prodAreaInt = nuclearInt .* cellArea;
% figure, histogram(prodAreaInt(prodAreaInt~=0),[],0);
% [~, cutoffValue] = cutFirstHistMode(prodAreaInt(prodAreaInt~=0));
% indxGood = find(prodAreaInt>=cutoffValue);
% indxBad = find(prodAreaInt<cutoffValue);
% 
% %get threshold based on intensity and area
% indxGood = find(nuclearInt>0);
% figure, histogram(cellArea(indxGood),[],0);
% [~, cutoffValue] = cutFirstHistMode(cellArea(indxGood));
% indxGoodGood = find(cellArea(indxGood)>cutoffValue);
% indxGood = indxGood(indxGoodGood);
% indxBad = setdiff(1:max(L(:)),indxGood);
% 
% %keep only good regions
% tf = ismember(L,indxBad);
% indxPixel = find(tf);
% L2 = L;
% L2(indxPixel) = -1;
% 
% %display results
% image3Color2(:,:,3) = (imageNfilter-min(imageNfilter(:)))/(max(imageNfilter(:))-min(imageNfilter(:)));
% image3Color2(:,:,1) = (L2-min(L2(:)))/(max(L2(:))-min(L2(:)));
% figure, imshow(image3Color2,[]);


% [~, cutoffValue] = cutFirstHistMode(imageC_NMS(imageC_NMS~=0));
% imageC_binary = double(imageC_NMS>=cutoffValue);
% imageC_mult = imageC_binary .* imageCfilter;

% imageCfilter = (imageCfilter - min(imageCfilter(:)))/...
%     (max(imageCfilter(:)) - min(imageCfilter(:)));
% level = graythresh(imageCfilter);
% imageC_mult = double(imageCfilter > level);

% imageC_mult = imageCfilter;

%subtract them from each other
% imageDiff = (max(imageNfilter(:))-min(imageNfilter(:)))*imageC_mult - ...
%     (max(imageC_mult(:))-min(imageC_mult(:)))*imageNfilter;
% imageDiff = (max(maskNuclei(:))-min(maskNuclei(:)))*imageC_mult - ...
%     (max(imageC_mult(:))-min(imageC_mult(:)))*maskNuclei;

% imageDiff = imageC_mult;
