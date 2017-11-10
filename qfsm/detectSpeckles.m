function [cands locMax]=detectSpeckles(I,noiseParam,specParam,varargin)

% detectSpeckles detects speckles in a filtered image
%
% Plots PSFs on the positions of the primary speckles (found by locmax operator) and substructs
% them from the filtered data. Appled again on the resulting image, the locmax-operator finds 
% new (secondary) speckles (intensity and distance significance tests applied)
%
%
% SYNOPSIS   cands=detectSpeckles(I,strg,counter,noiseParam,Speckles)
%
% INPUT      I          :  filtered image
%            noiseParam :  noise parameters for statistical speckle selection
%            specParam  :  (1) contains information for the hierarchical level 
%                          (2) minimal increase in (%) of new speckles
%                          (before stopping)
%
% OUTPUT     cands      :  augmented cands structure (see fsmPrepTestLocalMaxima.m)
%            locMax     :  image containing the significant local maxima 
%
% References:
% A. Ponti et al., Biophysical J., 84 336-3352, 2003.
% A. Ponti et al., Biophysical J., 89 2459-3469, 2005.

% Aaron Ponti, October 4th, 2002
% Sylvain Berlemont, Jan 2010
% Sebastien Besson, May 2011 (last modified Nov 2011)
% Adapted from fsmPrepMainSecondarySpeckles.m

ip=inputParser;
ip.addRequired('I',@isnumeric);
ip.addRequired('noiseParam',@isnumeric)
ip.addRequired('specParam',@isnumeric);
ip.addOptional('psfSigma',1,@isscalar);
ip.addParamValue('userROIbw',[],@isnumeric)
ip.parse(I,noiseParam,specParam,varargin{:});
userROIbw = ip.Results.userROIbw;
psfSigma = ip.Results.psfSigma;
IG=I;

% SB: What the hell is this constant!!!
SIG=1.88; % for the twice convolved image (or 1.77)

% local minima
Imin=locmin2d(IG,[3,3]);

% Add virtual points along the cell edge to make sure there will be a
% triangle in the delaunay triangulation for each local maxima (cands). If
% no cell mask is provided there is no additional point.

% reconstruct the mask from the original image
bwMask = I ~= 0;

if ~all(bwMask)
    % Dilate the cell mask
    bwMaskExt = imdilate(bwMask, strel('disk', 1));

    % Compute Distance transform
    bwDist = bwdist(1 - bwMaskExt);
    % Get first outline
    outline = contourc(double(bwDist), [0, 0]);
    X = [];
    Y = [];
    iChunk = 1;
    while iChunk < size(outline,2)
        n = outline(2,iChunk);
        X = [X, round(outline(1,iChunk+1:5:iChunk+n))];
        Y = [Y, round(outline(2,iChunk+1:5:iChunk+n))];
        iChunk = iChunk + n + 1;
    end
    % add these points to Imin
    ind = sub2ind(size(Imin), Y, X);
    Imin(ind) = -1000;
end

% intial (filtered) image
[cands,triMin,pMin] = fsmPrepConfirmSpeckles(IG,Imin,noiseParam,userROIbw);

[cands(:).speckleType] = deal(1);

Inew=IG;
candsS=cands;
HierLevel=2;

while HierLevel<=specParam(1) && ...
        length(candsS)>(specParam(2)*length(cands)) && ...
        any([candsS.status])
    
    Inew=fsmPrepSubstructMaxima(Inew,SIG,candsS);
    candsS = fsmPrepConfirmLoopSpeckles(Inew,noiseParam,triMin,pMin,IG,userROIbw);
    
    [candsS(:).speckleType] = deal(HierLevel);
    
    candsS=fsmPrepCheckDistance(candsS,cands,psfSigma);
    
    HierLevel=HierLevel+1;
    
    if ~isempty(candsS)
        cands=cat(1,cands,candsS);
    end
end

% remove repetitions because of secondary speckles apearing on the same
% positions as primary (because of floating background)
cands=removeSpeckleRepetitions(cands);

% obtain updated IM from candstTot
locMax=zeros(size(IG));

validCands = [cands(:).status] == 1;
if any(validCands)
    validLmax = vertcat(cands(validCands).Lmax);
    validILmax = [cands(validCands).ILmax];
    validIdx = sub2ind(size(IG), validLmax(:,1), validLmax(:,2));
    locMax(validIdx) = validILmax; 
end

function candsTot=removeSpeckleRepetitions(candsTot)
% removeSpeckleRepetitions checks for repetitions amongst the lists of the primary, 
% secondary and tertiary speckles and removes them from the cands structure
%
% Sebastien Besson, June 2011
% Adapted from fsmPrepCheckInfo

candsTotPos = vertcat(candsTot.Lmax);
if ~isempty(candsTotPos),
    % Construct kd-tree
    idx=KDTreeBallQuery(candsTotPos,candsTotPos,0);
else
    idx = {};
end
% Find replicates and return if no duplicate found
isReplicate = cellfun(@(x) length(x)>1,idx);
if isempty(find(isReplicate,1)), return; end

% Check that no two candidates are significant
checkCands = find(cellfun(@(x) sum([candsTot(x).status])>1, idx(isReplicate)),1);
if ~isempty(checkCands),
    error('unusual speckles repetition: both speckles (at the same position) are significant');
end

% Identify insignificant maxima and remove them
ind =([candsTot.status]==0 & [candsTot.speckleType]>1 & isReplicate');
candsTot(ind)=[];

function newCands=fsmPrepCheckDistance(newCands,oldCands,psfSigma)

% fsmPrepCheckDistance uses a look-up-table to verify if a new (secondary or tertiary) speckle are not closer than 
% theoretically possible to a primary (for tertiary also secondary) one; if a candidate is closer - it is discarded
%
% SYNOPSIS    newCands=fsmPrepCheckDistance(newCands,oldCands)
%
% INPUT      newCands      :   cands for the new/secondary speckles
%            oldCands      :   cands for the old/primary speckles
%            psfSigma      :   standard deciation of the Gaussian psf
%
% OUTPUT     newCands      :   updated cands for the new/secondary speckles
%
%
% DEPENDENCES   fsmPrepCheckDistance uses { createDistanceMatrix }
%               fsmPrepCheckDistance is used by { fsmPrepMainSecondarySpeckles }
%
% Alexandre Matov, November 7th, 2002
% Modified by Sylvain Berlemont, 2010
% Modified by Sebastien Besson, Oct 2011

% Retrieve only statistically significant speckles
validIdx1 = [oldCands(:).status] == 1;
validIdx2 = [newCands(:).status] == 1;
if ~(any(validIdx1) && any(validIdx2)),return; end

% Get positions and intensities of i-th and i+1th order speckles
pos1 = vertcat(oldCands(validIdx1).Lmax);
pos2 = vertcat(newCands(validIdx2).Lmax);

deltaI1 = [oldCands(validIdx1).deltaI];
deltaI2 = [newCands(validIdx2).deltaI];

% Old distance criterion from fsmCenter
% D=createSparseDistanceMatrix(pos2,pos1,5.21);
% toBeRemoved=false(size(D,1),1); % initialization
% 
% for i=1:size(D,1)
%     rowI = D(i,:);
%     minDistRow=find(rowI);
%     
%     if ~isempty(minDistRow)
%         % SB: this criteria doesn't make sense to me:
%         RelInt = deltaI1(minDistRow) / deltaI2(i);
%         r=5.2-3.5.*exp(.5-RelInt(:));
%         d = nonzeros(rowI) + eps;
%             
%         toBeRemoved(i) = any(d <= r);
%     end
% end

% Classify secondary speckles by their distance to the primary speckles
rmax = 4*psfSigma;
idx1= KDTreeBallQuery(pos1,pos2,2*psfSigma);
[idx2 D2]= KDTreeBallQuery(pos1,pos2,rmax);

% Intialize logical array for removing speckles
r=false(sum(validIdx2),1);

% Remove all points closer than 2*psfSigma to any existing speckle
r(~cellfun(@isempty,idx1))=true;

% Apply generalized Sparrow criterion to all speckles which distance is
% between 2 psfSigma and 4 psfSigma of an existing speckle
% See A. Ponti et al. 2005 (supplementary material)
w = @(r) sqrt(r.^2-4*psfSigma^2);
Acrit = @(r) (r-w(r))./(r+w(r)).*exp(w(r).*r/(2*psfSigma.^2));
A = @(x) deltaI1(idx2{x}) / deltaI2(x);
applyGSC = @(x) all(A(x)<=Acrit(D2{x})');

for i=find(~cellfun(@isempty,idx2) & cellfun(@isempty,idx1))'
    r(i)=~applyGSC(i);
end

% Remove invalid speckles
if any(r)
    ind2 = find(validIdx2);
    newCands(ind2(r)) = [];
end

 function [cands,triMin,pMin]=fsmPrepConfirmSpeckles(IG,Imin,noiseParam,userROIbw)
% fsmPrepConfirmSpeckles uses statistical tests to confirm the significance
% of detected speckles
%
% SYNOPSIS   [cands,triMin,pMin]=fsmPrepConfirmSpeckles(IG,Imin,noiseParam)
%
% INPUT      IG         :  filtered image
%            Imin       :  local minima
%            noiseParam :  noise parameters for statistical speckle selection
%            userROIbw  :  (optional) User-defined black-and-white mask to select speckles
%
% OUTPUT     cands      :  cands structure (see fsmPrepTestLocalMaxima.m)
%            triMin     :  set of Delaunay triangles
%            pMin       :  attached dSet (vertex coordinates) to lMin
%
%
%
% DEPENDENCES   fsmPrepConfirmSpeckles uses { fsmPrepBkgEstimationDelaunay, fsmPrepTestLocalMaxima, locmax2d }
%               fsmPrepConfirmSpeckles is used by { fsmPrepMainSecondarySpeckles }
%
% Alexandre Matov, November 7th, 2002
% Sylvain Berlemont, Jan 2010

% Find the local maxima  
Imax=locmax2d(IG,[5,5]);

if nargin == 4 && ~isempty(userROIbw) 
    Imax=Imax.*userROIbw;
end

% Finds 3 loc min around each loc max
[cands,triMin,pMin]=fsmPrepBkgEstimationDelaunay(Imax,Imin);

% analyze speckles - validate, locmax, locmin...
cands = fsmPrepTestLocalMaxima(IG,cands,noiseParam,IG);  


function cands=fsmPrepConfirmLoopSpeckles(Inew,noiseParam,triMin,pMin,IG,userROIbw)

% fsmPrepConfirmLoopSpeckles uses statistical tests to confirm the significance of detected speckles
% of higher than one hierarchical level (in the main loop)
%
% SYNOPSIS   cands=fsmPrepConfirmLoopSpeckles(Inew,noiseParam,triMin,pMin,IG)
%
% INPUT      Inew       :  substracted image
%            noiseParam :  noise parameters for statistical speckle selection
%            triMin     :  set of Delaunay triangles
%            pMin       :  attached dSet (vertex coordinates) to lMin
%            IG         :  original (filtered) image
%
% OUTPUT     cands      :  cands structure (see fsmPrepTestLocalMaxima.m)
%
%
%
% DEPENDENCES   fsmPrepConfirmLoopSpeckles uses { locmax2d, tsearch,fsmPrepTestLocalMaxima }
%               fsmPrepConfirmLoopSpeckles is used by { fsmPrepMainSecondarySpeckles }
%
% Alexandre Matov, April 2nd, 2003
% Sylvain Berlemont, Jan 2010

% find the local maxima  
Imax=locmax2d(Inew,[5,5]);

if nargin == 6 && ~isempty(userROIbw)
    % Mask Imax
    Imax=Imax.*userROIbw;
end

% find the coordinates/positions of the initial local maxima before
% significance test (for comparision)
[yi,xi]=find(ne(Imax,0));
pMax=[yi,xi];

% Assign local maxima to local minimum triangles
% triangles=tsearch(pMin(:,1),pMin(:,2),triMin,pMax(:,1),pMax(:,2));
triangles = triMin.pointLocation(pMax);

% Store information into cands structure
validTri = 1:numel(triangles);
validTri = validTri(~isnan(triangles));
n = ones(numel(validTri), 1);

cands = struct(...
    'Lmax', mat2cell(pMax(validTri,:), n), ...
    'Bkg1', mat2cell(pMin(triMin(triangles(validTri),1),:), n),...
    'Bkg2', mat2cell(pMin(triMin(triangles(validTri),2),:), n),...
    'Bkg3', mat2cell(pMin(triMin(triangles(validTri),3),:), n));

% analyze speckles - validate, locmax, locmin...
cands = fsmPrepTestLocalMaxima(Inew,cands,noiseParam,IG);  


function [Inew,Imaxima]=fsmPrepSubstructMaxima(IG,SIG,cands)

% fsmPrepSubstructMaxima takes the filtered image IG and substracts Gaussian
% kernels at the positions of the local maxima
%
%
% SYNOPSIS   Inew=fsmPrepSubstructMaxima(IG,Imax,SIG,cands)
%
% INPUT      IG         :  filtered image
%            SIG        :  sigma of the GKs
%            cands      :  cands strucure
%
% OUTPUT     Inew       :  the resulting image after substaction
%            Imaxima    :  the synthetic image to be substracted from the
%                          real
%
% Alexandre Matov, November 7th, 2002
% Modified by Sylvain Berlemont, 2010

Imax=zeros(size(IG)); % prepearing Imax for substraction

validCands = ([cands(:).status] == 1) & ([cands(:).deltaI] > 0);
if any(validCands)
    validLmax = vertcat(cands(validCands).Lmax);
    validDeltaI = [cands(validCands).deltaI];
    validIdx = sub2ind(size(IG), validLmax(:,1), validLmax(:,2));
    Imax(validIdx) = validDeltaI;
end

% the masked with a GK local maxima points with intensity delta.I for the
% raw data image.
Imaxima=filterGauss2D(Imax,SIG);

% substruction
Inew=IG-Imaxima;

% In case Imaxima has value outside the cell footprint (IG is masked) due
% to the filtering step, Inew may have negative value outside the cell
% footprint. Make sure the cell mask is also applied on Inew.
Inew(IG == 0) = 0;

function [cands,triMin,pMin]=fsmPrepBkgEstimationDelaunay(lMax,lMin)
% fsmPrepBkgEstimationDelaunay uses Delaunay triangulation to assign 3 local minima to every local maximum
%
% SYNOPSIS   [cands,triMin,pMin]=fsmPrepBkgEstimationDelaunay(lMax,lMin)
%
% INPUT      lMax      :   local max map (the output of the locMax2D function)
%            lMin      :   local min max (the output of the locMin2D function)
%
% OUTPUT     triMin    :   Delaunay triangulation
%            pMin      :   attached dSet (vertex coordinates) to lMin
%            cands     :   structure containing statistical information for each local maximum
%                          (see help for detail) - fsmPrepBkgEstimationDelaunay only stores 
%                          local maximum and local minimum coordinates
%
% DEPENDENCES
%
% Aaron Ponti, October 23th, 2002
% Sylvain Berlemont, Jan 2010
% Sebastien Besson, May 2011

% Collect lMax coordinates into matrics
[y x]=find(lMax);
pMax=[y x];

if isempty(pMax), return; end

% Collect lMin coordinates into matrices
[y x]=find(lMin);
pMin=[y x];

% Delaunay triangulation
triMin=DelaunayTri(pMin);
triangles = triMin.pointLocation(pMax);

% Store information
validTri = 1:numel(triangles);
validTri = validTri(~isnan(triangles));
n = ones(numel(validTri), 1);

cands = struct(...
    'Lmax', mat2cell(pMax(validTri,:), n), ...
    'Bkg1', mat2cell(pMin(triMin(triangles(validTri),1),:), n),...
    'Bkg2', mat2cell(pMin(triMin(triangles(validTri),2),:), n),...
    'Bkg3', mat2cell(pMin(triMin(triangles(validTri),3),:), n));

function cands = fsmPrepTestLocalMaxima(img,cands,parameters,imgO)
% fspPrepTestLocalMaxima selects statistically significant local maxima using loc max and background info from analyzeSpeckles
%
% SYNOPSIS   cands=fsmPrepTestLocalMaxima(img,cands,parameters)
%
% INPUT      img        : image matrix
%            cands      : Imax and Imin for every speckles obtained 
%                         from analyzeSpeckles (Delaunay triangulation)
%            parameters : [k sigmaD PoissonNoise I0]
%                         k : gives the confidence interval Xavg+/-k*sigma 
%                         for normally distributed data
%                         sigmaD, PoissonNoise, I0 : from the noise
%                         calibration model
%                         sigma=sqrt(sigmaD^2+PoissonNoise*(I-I0))
%                         I0 : mean intensity of the background images
%            imgO       : the original (filtered) image 
%                              
% OUTPUT     cands      : cands (input) augmented by additional (statistical) information 
%
%                         cands.Lmax       : Local maximum position - [y x]
%                              .Bkg1       : First local minimum position - [y x]
%                              .Bkg2       : Second local minimum position - [y x]
%                              .Bkg3       : Third local minimum position - [y x]
%                              .ILmax      : Local maximum intensity
%                              .IBkg       : Mean background intensity
%                              .deltaI     : Intensity difference: ILmax-IBkg
%                              .deltaICrit : Critical intensity difference as calculated with the noise model
%                              .sigmaLmax  : Error on local maximum intensity
%                              .sigmaBkg   : Error on background intensity 
%                              .status     : Significance of the local maximum: 1, speckle; 0, weak local maximum
%

% Aaron Ponti, October 4th, 2002
% Sylvain Berlemont, Jan 2010
% Sebastien Besson, Aug 2011

% Parameters
k=parameters(1);
sigmaD=parameters(2);
PoissonNoise=parameters(3);
I0=parameters(4);

Lmax = vertcat(cands(:).Lmax);
Bkg1 = vertcat(cands(:).Bkg1);
Bkg2 = vertcat(cands(:).Bkg2);
Bkg3 = vertcat(cands(:).Bkg3);

indLmax = sub2ind(size(img), Lmax(:, 1), Lmax(:, 2));
indBkg1 = sub2ind(size(img), Bkg1(:, 1), Bkg1(:, 2));
indBkg2 = sub2ind(size(img), Bkg2(:, 1), Bkg2(:, 2));
indBkg3 = sub2ind(size(img), Bkg3(:, 1), Bkg3(:, 2));

% Calculate difference Imax-Imin for every speckle
for c1=1:length(cands)
	
    % Read values for Imax and Imin from the image at the coordinates
    % specified by info 
    Imax = img(indLmax(c1));
    
    % imgO is the original filtered image
    Imin = mean(nonzeros([imgO(indBkg1(c1)), imgO(indBkg2(c1)), imgO(indBkg3(c1))]));
    
    [Imin,deltaI,k,sigmaDiff,sigmaMax,sigmaMin,status]=...
        testSpeckleSignificance(Imax,Imin,k,sigmaD,PoissonNoise,I0);
	
	% Complete cands
    cands(c1).ILmax=Imax;             % Local maximum intensity
    cands(c1).IBkg=Imin;              % Mean background intensity
    cands(c1).deltaI=deltaI;          % Intensity difference: ILmax-IBkg
    cands(c1).deltaICrit=k*sigmaDiff; % Critical intensity difference as calculated with the noise model
    cands(c1).sigmaLmax=sigmaMax;     % Error on local maximum intensity
    cands(c1).sigmaBkg=sigmaMin;      % Error on background intensity 
    cands(c1).status=status;          % Significance of the local maximum: 1, speckle; 0, weak local maximum
end

% Filter out cands with no local background
cands(isnan([cands.IBkg]))=[];