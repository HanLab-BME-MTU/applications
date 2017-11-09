function [cT, cNull, cCritical] = colocalMeasurePt2Blob(SegmentationCoords,detectionCoords, threshold, maskingFile,varargin)           
%COLOCALMEASUREPT2Blob measures the colocalization between segmentation and detection point
%
%Synopsis: [cT, cNull, cCritical] = colocalMeasurePtBlob(SegemntationCoords,...
%                      detectionObs, threshold,maskingFile,labels,varargin)
%
%   Function analyzes the distances of nearest neighbors between detection 
%   point and blob segmentation and uses a distance threshold to ouput the 
%   probability of finding points in the observed channel that colocalize 
%   with the blob segementation, the associated interaction potential 
%   (see source below) and the null hypothesis values for both of these.
%Input:
% SegmentationCoords: Coordinates of the segmentations.
% 
%
% detectionCoords: Detected positions of molecules in the observed channel
% 
%
% Threshold: distance threshold for possible interacting objects
%
% Masking File: This should either be a binary mask, or if no mask is given 
%               then the size of the image should be input to create a mask of all ones
%
% varargin (labels): labels from segmentations.
%                    Optional. Only needed for get detections inside blobs
%
% Output:
% cNull (Null Hypothesis): probabilty of finding a distance <= threshold to
%   a detection in random detections under the null hypothesis (No interaction)
% 
% cT (Colocalization Measure): probabilty of finding a detection in detectionObs
%   a distance <= threshold to a random detection
%
% cCritical: minimum probability needed to consider cT significant
%
% Estimator methodology adapted from "Beyond co-localization: inferring
% spatial interactions between sub-cellular structures from microscopy
% images" - Helmuth et al. BMC Bioinformatics 2010 
% 
% Anthony Vega 09/2014
%Jesus Vega July 2017


%% Input

%Observed detections
  xdetec = detectionCoords.xCoord(:,1);
  ydetec = detectionCoords.yCoord(:,1);
 
 Detections = [ydetec xdetec];

%Find detections inside the mask
if isempty(maskingFile)
    error('A binary mask or size of image should be inputed')

elseif ~ismatrix(maskingFile)
    maskingFile = ones(maskingFile);
end

[row, col] = find(maskingFile); 
    maskList = [row, col];
    lia1 = ismember(round(Detections),maskList,'rows');
    
%Multipling lia (binary vector) by Detections will replace coord outside
    %boundary with zero, last line removes all zeros from vector
    Detections(:,1) = lia1.*Detections(:,1);
    Detections(:,2) = lia1.*Detections(:,2);
    Detections( ~any(Detections,2), : ) = [];
    

%Segementations
%Find segmentations inside the mask
    lia1 = ismember(round(SegmentationCoords),maskList,'rows');
    
%Multipling lia (binary vector) by Segmentations will replace coord outside
    %boundary with zero, last line removes all zeros from vector
    SegmentationCoords(:,1) = lia1.*SegmentationCoords(:,1);
    SegmentationCoords(:,2) = lia1.*SegmentationCoords(:,2);
    SegmentationCoords( ~any(SegmentationCoords,2), : ) = [];


%Get nearest neighbor of observed detections
[~, NNDist] = knnsearch(SegmentationCoords,Detections,'k',1);

%% Test NN under threshold
 %Probability density of observed detections
 [f,xi] = ksdensity(NNDist);    
  
%Determine density of distances below threshold to find cNull
%commented lines take distances from the periphery
 %testmin = find(xi >= 1);
 %testmax = find(xi <= threshold);
 %test = testmin(ismember(testmin,testmax));
 test = find(xi<=threshold);
 cT = trapz(xi(1:max(test)),f(1:max(test))); 

%% Detections Inside Blobs

if ~isempty(varargin)
%Find detections inside the blobs
BigBlobLoc = ismember(SegmentationCoords,round(Detections),'rows');

%Get the coordinates of the points where is a detections inside the blob
NNCoord = zeros(size(SegmentationCoords));
NNCoord(:,1) = BigBlobLoc.*SegmentationCoords(:,1);
NNCoord(:,2) = BigBlobLoc.*SegmentationCoords(:,2);
NNCoord(~any(NNCoord,2),:) = [];

%Convert coordinate to Idxs                                          
NNidx = sub2ind(size(maskingFile),NNCoord(:,1),NNCoord(:,2));

maskforDetec = zeros(512);
maskforDetec(NNidx) = 1;

%Get the position of detections inside the blob 
findDetec = find(maskforDetec);

%Get the position of the blobs             
findLabels = find(varargin{1});

%Get the label of the blob where is a detection
DetecInLabel = ismember(findLabels,findDetec);

%Get the number of the label
LabelNum  = varargin{1}(findLabels(DetecInLabel));

%Pick the label that has a detection
final = ismember(varargin{1},LabelNum);

%plot detections inside the blobs
justborder = bwmorph(final,'remove');
figure('Name','Blobs with Detections Inside'), imshowpair(justborder,maskforDetec);
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %Intensities and eccentricities of segmentations with detections inside
% SegIntensity = regionprops(final,imageNorm,'PixelValues');
% 
% %SegEccentricity = regionprops(final,'Eccentricity');
% 
% %Display the number of segmentations with detections inside
% NumOfBlobsWithDetec  = numel(SegIntensity);
% %disp(['Number of Segementations with Detections Inside: ' num2str(NumOfBlobsWithDetec)]);  
% 
% %Calculate inegrated intensity for each blob with detections inside
% Integrated_Int = zeros(NumOfBlobsWithDetec,1);
% 
% for k = 1:NumOfBlobsWithDetec
%     
%     Integrated_Int(k) = sum(SegIntensity(k).PixelValues);
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%PLOTTING

% %plot all blobs
% Segidx = sub2ind(size(maskingFile),SegmentationCoords(:,1),SegmentationCoords(:,2));
% mask = zeros(512);
% mask(Segidx) = 1;
% figure('Name','All Blobs'); imshow(mask,[]);
%
% 
% %plot the areas of the blobs with detections inside
%   BlobSize = regionprops(final,'Area');
%   BlobSize = cell2mat(struct2cell(BlobSize))';
%   
%  figure, histogram(BlobSize,20:20:(round(max(BlobSize),-2)+60));

   
%% Test against null hypothesis

[~, RandDist] = knnsearch(SegmentationCoords,maskList,'K',1);

%Probability density of 
[f,xi] = ksdensity(RandDist);
    
%Determine density of distances below threshold to find cNull
%commented lines take distances from the periphery
% testmin = find(xi >= 0);
% testmax = find(xi <= threshold);
%test = testmin(ismember(testmin,testmax));
test = find(xi<=threshold);
cNull(1,1) = trapz(xi(1:max(test)),f(1:max(test)));

% Estimate critical parameters based on cNull and interaction size
cCritical = (binoinv(0.95,length(NNDist),cNull(1,1)))/length(NNDist);
end

