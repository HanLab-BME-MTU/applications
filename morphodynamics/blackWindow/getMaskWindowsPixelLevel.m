function winMat = getMaskWindowsPixelLevel(mask,perpSize,paraSize,varargin)
%Wrting the windowing the way I should have years ago... God damnit I can't
%believe I didn't do it this way before ARRRRRGGGGG!!!

%TEMP!
%FINISH DOCUMENTATION

%% ---------------------- Input ------------------------------ %%

ip = inputParser;
ip.FunctionName = mfilename;
ip.KeepUnmatched = true;%Keep unmatched parameters so we can warn the user
ip.addRequired('mask',@(x)(islogical(x) && ndims(x) == 2));
ip.addRequired('perpSize',@(x)(numel(x) == 1 && isfinite(x) && x >= 1));
ip.addOptional('paraSize',[],@(x)(numel(x) == 1 && isfinite(x)));
ip.addParamValue('StartPoint',[],@(x)(ndims(x) == 2 && size(x,2)  == 2))
%ip.addParamValue('StartContour',1,@(x)(numel(x) == 1 && isposint(x)));

%TEMP - is this what we want to do? This creates a field in p which has the
%mask and perpSize in it also, even though we have variables for those!!
%Wastes some small amount of memory...
ip.parse(mask,perpSize,varargin{:});
p = ip.Results;

if ~isempty(fieldnames(ip.Unmatched))
    warning(['MORPHODYNAMICS:' mfilename ':unrecognizedParameter'],'Unrecognized parameters were input! These parameters will be ignored!');
end

showPlots = false;

%ADDITIONAL INPUT CHECKING!!??!
%need to check image size and switch to bwdist_old if it's too big...

%% ------------------------ Init ----------------------------- %%

%Get distance transform and label matrix
[distX,distL] = bwdist(~mask);
maxDist = max(distX(:));
[M,N] = size(mask);

% ---------- Distance Isocontour Values -------- %
%Determine the distance values bounding each "band" which will divide the
%mask in the direction perpindicular to the edge

distXvals = 0:perpSize:maxDist;
if distXvals(end) ~= maxDist
    distXvals = [distXvals maxDist];
end
nBands = numel(distXvals)-1;

% -------- Start Contour Creation ----------- %
%Gets the pixels of the "Start Contour" which is the contour where the
%windows have the size specified by paraSize

%Use contourc to get the pixel locations because it is much faster than
%bwboundaries...
startContour = separateContours(contourc(double(mask),[1 1]));
%TEMP - do clean-up but we 


%And get the linear indices corresponding to these pixels
startContour = startContour{1}';
iStartContour = sub2ind(size(mask),startContour(:,2),startContour(:,1));

%TEMP - faster to do this logically since we only have 1 and sqrt(2)
%values?
dAlong = [0 cumsum(sqrt(sum(diff(startContour(:,1),1,1) .^2 + diff(startContour(:,2),1,1) .^2,2)))'];
dVals = 0:paraSize:dAlong(end);
nSlices = numel(dVals)-1;
[~,iBestDist] = arrayfun(@(x)(min(abs(dAlong - dVals(x)))),1:nSlices+1);

if showPlots
   figure
   imagesc(distL),axis image,colormap gray,hold on
   plotDirection(startContour);
   plot(startContour(iBestDist,1),startContour(iBestDist,2),'rx');
    
end

%TEMP - adjust class depending on nWinMax...
winMat = zeros(size(mask),'uint16');

for j = 2:nSlices
    
    %Get indices of border pixels for this slice
    iCurrBord = iStartContour(iBestDist(j-1):iBestDist(j));
    
    %And indices of closest background pixels corresponding to these
    iCurrBack = distL(iCurrBord);
    
    %Now find all pixels which share these labels
    mSlice = reshape(any(bsxfun(@eq,distL(:),iCurrBack'),2),size(mask));
    
    winMat(mSlice) = j-1;
    
end

winMat(~mask) = 0;

dPerpVals = 0:perpSize:max(distX(:));
if dPerpVals(end) ~= max(distX(:))
    dPerpVals = [dPerpVals max(distX(:))];
end

nBand = numel(dPerpVals)-1;

%The first band keeps their original indices
for j = 3:nBand+1
   
    currBand = (distX>dPerpVals(j-1) & distX <= dPerpVals(j));
    winMat(currBand) = winMat(currBand) + nSlices*(j-2);
    
    
    
end






