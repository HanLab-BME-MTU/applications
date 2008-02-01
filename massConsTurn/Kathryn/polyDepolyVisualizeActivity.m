function polyDepolyVisualizeActivity(cmDir,mapDir,edgePix,gamma,cMapLength,batchGlobalMinMax)
% POLYDEPOLYVISUALIZEACTIVITY uses isomorphic color maps to show poly/depoly activity
%
% DESCRIPTION: uses isomorphic color maps to show poly/depoly activity
%
% SYNOPSIS: polyDepolyVisualizeActivity(cmDir,mapDir,edgePix,gamma,cMapLength,batchGlobalMinMax)
%
% INPUT: 
%    cmDir             : path to analysis\edge\cell_mask
%    mapDir            : path to ONE LEVEL UP FROM analysis\turn\...\mapMats
%    edgePix           : .mat saved in analysis\turn by polyDepolyMap
%    gamma             : value for gamma correction (<1 boosts low signals)
%                        (default = 1: no change)
%    cMapLength        : color map length for red/green (default is 64).
%                        must be even
%    batchGlobalMinMax : [gmin gmax] this optional parameter is useful if
%                        you want to normalize a bunch of movies with the 
%                        same min/max values. if empty or 0, the function
%                        uses the values specific to this movie.
%                     
% OUTPUT: 
%    a new directory is created at the level of \turn\...\mapTifs where
%    output tifs are stored.  several parameters are saved in the .mat file
%    visualizationParams
%
% Note on batchGlobalMinMax: if running a batch job and you use the batch
% values, and if one movie has a very broad range compared to the others,
% all the other movies will tend to gray.  this is because that movie will
% have optimized red/green color and the others will be closer to zero.
%
% It will therefore be best when running as a batch to set a good range
% from an empirical estimate of the activity histograms made during time
% averaging, even if this cuts off the broad range values to the
% boundaries.
%
% MATLAB VERSION (originally written on): 7.2.0.232 (R2006a) Windows XP
%
% USERNAME: kathomps
% DATE: 12-Apr-2006
%

if nargin<1 || ~isdir(cmDir)
    error('polyDepolyVisualizeActivity: cell mask directory does not exist')
end

if nargin<2 || ~isdir([mapDir filesep 'mapMats'])
    error('polyDepolyVisualizeActivity: mapDir should contain mapMats directory')
end

if nargin<3 || ~iscell(edgePix)
    error('polyDepolyVisualizeActivity: edgePix required (saved under turn directory)')
end

if nargin<4 || isempty(gamma)
    gamma=1; % default - no correction
    disp('polyDepolyVisualizeActivity: gamma=1 (no gamma correction)')
elseif gamma<=0
    error('polyDepolyVisualizeActivity: gamma must be greater than 0')
end

if nargin<5 || isempty(cMapLength)
    cMapLength=64; % default
    disp('polyDepolyVisualizeActivity: default cMapLength=64')
elseif ~isEven(cMapLength)
    error('polyDepolyVisualizeActivity: cMapLength should be even')
end

if nargin<6 || isempty(batchGlobalMinMax) || sum(abs(batchGlobalMinMax))==0 % default - will use actual max/min from this movie
	batchGlobalMinMax=0;
    disp('polyDepolyVisualizeActivity: will use values from this movie for batchGlobalMinMax');
elseif length(batchGlobalMinMax)~=2
    error('polyDepolyVisualizeActivity: batchGlobalMinMax should be [] or [userMin userMax]')
end


% get list of activity maps and count them
[listOfMaps] = searchFiles('polyDepoly','tif',mapDir,1);
nFrames=size(listOfMaps,1);

s=length(num2str(nFrames));
strg=sprintf('%%.%dd',s);

% get list of corresponding cell masks
[listOfCellMasks] = searchFiles('.tif',[],cmDir,0);

% get min/max over whole image series if not using global values from batch
if batchGlobalMinMax==0
    minValue=0;
    maxValue=0;
    for i=1:nFrames
        indxStr=sprintf(strg,i);
        iMap=load([char(listOfMaps(i,2)) filesep char(listOfMaps(i,1))]);
        im=eval(['iMap.' char(fieldnames(iMap))]);

        if i==1
            [imL,imW]=size(im);
        end

        minValue=min(minValue,nanmin(im(:)));
        maxValue=max(maxValue,nanmax(im(:)));
    end
else
    minValue=batchGlobalMinMax(1);
    maxValue=batchGlobalMinMax(2);
end

% if gamma correction, apply it to the min/max values so scaling is right
if gamma~=1
    sign1=1;
    sign1(minValue<0)=-1;
    minValue=sign1.*(abs(minValue).^gamma);
    
    sign2=1;
    sign2(maxValue<0)=-1;
    maxValue=sign2.*(abs(maxValue).^gamma);
end

% m is the upper bound for normalization to color map range
% here we take the minimum abs value and will make values < -m equal to -m
% and values > m equal to m

% this means if running a batch job and one movie has a very broad range
% compared to the others, all the other movies will tend to gray, b/c that
% movie will have optimized red/green color and the others will be closer
% to zero.

% it will therefore be best when running as a batch to set a good range
% from an empirical estimate of the activity histograms made during time
% averaging, even if this cuts off the broad range values to the
% boundaries.

m=min(abs([minValue maxValue]));

% make output directory
outDir=[mapDir filesep 'mapTifs_' num2str(m) '_' num2str(gamma)];
if ~isdir(outDir)
    mkdir(outDir);
else
    delete([outDir filesep '*tif'])
end


% create the color map with black at the end
% this function makes a perceptual color map
colorMap=isomorphicColormap('g/r',cMapLength);
colorMap=[colorMap; [0 0 0]];


for i=1:nFrames
    % load the map (im) and cell mask (cellMask)
    indxStr=sprintf(strg,i);
    iMap=load([char(listOfMaps(i,2)) filesep char(listOfMaps(i,1))]);
    im=eval(['iMap.' char(fieldnames(iMap))]);

    if i==1
        [imL,imW]=size(im);
    end

    % make mask into RGB-size array
    cMask=double(imread([char(listOfCellMasks(i,2)) filesep char(listOfCellMasks(i,1))]));
    cellMask=repmat(cMask,[1 1 3]); 

    % gamma-correct pixel values here (boost low values if gamma < 1)
    sign=zeros(size(im));
    sign(im<0)=-1;
    sign(im>=0)=1;
    im=sign.*(abs(im).^gamma);

    % rescale the image
    im(im<-m)=-m;
    im(im>m)=m;
    im=im+m;                    % range is 0-2m
    im=im./(2*m);               % range is 0-1
    im=im.*(cMapLength-1);      % range is 0-63 (if cMapLength is 64)
    im=round(im);               % round to integers
    im=im+1;                    % range is 1-64 (-32 to 32)
    im(isnan(im))=cMapLength+1;

    % initialize RGB arrays
    R=zeros(imL,imW);
    G=zeros(imL,imW);
    B=zeros(imL,imW);

    % fill RGB arrays with image values
    R(:)=colorMap(im(:),1);
    G(:)=colorMap(im(:),2);
    B(:)=colorMap(im(:),3);

    % make the cell edge white
    R(edgePix{i})=1;
    G(edgePix{i})=1;
    B(edgePix{i})=1;

    % make the final RGB image
    RGmap=zeros(imL,imW,3);
    RGmap(:,:,1)=R;
    RGmap(:,:,2)=G;
    RGmap(:,:,3)=B;

    % apply the cell mask
    RGmap=RGmap.*cellMask;

    % write the image
    imwrite(RGmap,[outDir filesep 'polyDepoly' indxStr '.tif']);
    % disp(['Saving red-green tiffs: frame ' num2str(i)])

end
save([outDir filesep 'visualizationParams'],'minValue','maxValue','m','gamma','cMapLength','batchGlobalMinMax');
