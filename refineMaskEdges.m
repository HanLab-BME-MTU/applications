function refinedMask = refineMaskEdges(maskIn,imageIn,maxAdjust,maxGap,preGrow)


% refinedMask = refineMaskEdges( mask , image , maxEdgeAdjust, maxEdgeGap,preGrow)
%
%This function uses edge detection to refine the edges of a mask.
%It assumes that the input mask maskIn OVERSHOOTS the cell boundary in
%all areas - that is the cell/object in the image is completely contained
%within the mask. If you aren't sure this is the case, input a non-zero
%integer for preGrow - this will grow the mask a little before refining it.
% the parameter maxEdgeAdjust is the maximum distance in pixels that the edge
% will be adjusted.
% maxEdgeGap is the the closure radius used to close the edges. This is
% effectively the largest gap that may occur in the edges in imageIn

%Hunter Elliott 3/2009


if nargin < 5 || isempty(preGrow)
    preGrow = 0;
elseif abs(round(preGrow)) ~= preGrow
    error('The mask growth radius, preGrow, must be a positive integer!')
end

if nargin < 4 || isempty(maxGap)
    maxGap = 5;
end

if nargin < 3 || isempty(maxAdjust)
    maxAdjust = 5;    
end

if ~isequal(size(maskIn),size(imageIn));
    error('The input mask and image must be the same size!!!')
end


showPlots = 0;

%Disable rounding warning
warning('off','MATLAB:intConvertOverflow');

%% ------ Parameters ----- %%


tooSmall = 30; %If a mask fragment is below this size in pixels it is thrown out.
threshScale = .6; %Fraction by which to adjust edge threshold on second round.
sigFilter = 1.0; %Sigma of filter used in canny edge detection.

%closeRad = max(round(min(size(imageIn)) / 5), 10); %Radius of closure op for edges

%% ----- Edge Detection ------ %%

%Run initial edge detection
[edges,autoThresh] = edge(imageIn,'canny',[],sigFilter);


%Run a second round with lower threshold
[edges,eThres] = edge(imageIn,'canny', autoThresh * threshScale,sigFilter);


%% --- Edge Refinement ---- %%

%Clean up the mask a little
maskIn = bwareaopen(maskIn,tooSmall);

%Grow the mask, if requested
if preGrow > 0 
    seGrow = strel('disk',preGrow);
    maskIn = imdilate(maskIn,seGrow);
end

%Get the distance transform of the mask and its inverse
distX = bwdist(~maskIn);

%Get rid of edges that are further than maxAdjust away from the mask edge,
%or which are outside of the edge.
edges = edges & distX <= maxAdjust & distX > 0;

%Add an inner border at maxAdjust inwards from the mask
edges = edges | bwperim(distX >= maxAdjust);

%Close these edges
seClose = strel('disk',maxGap);
closedEdges = imclose(edges,seClose);

if showPlots
    figure
    imagesc(imageIn)
    axis image, colormap gray,axis off
    hold on
    spy(bwperim(maskIn),'r',10)
    spy(closedEdges,'b',10)
    spy(edges,'y',10)
   
end

%Fill the center and any other holes in these edges to get the mask
refinedMask = imfill(closedEdges,'holes');

% %Keep only the largest object
% [labelCut,nCut] = bwlabeln(closedEdges,4);
% areas = zeros(1,nCut);
% for j = 1:nCut
%     
%     areas(j) = nnz(labelCut == j  & maskIn);
%     
% end

% %Find largest area
% [maxArea,iMax] = max(areas);
% 
% refinedMask = labelCut == iMax;


if showPlots
    
    figure
    imagesc(imageIn)
    axis image, colormap gray,axis off
    hold on
    spy(bwperim(maskIn),'r',15)    
    spy(bwperim(refinedMask),'b',10)
    %spy(edges,'y',5)
end
    


