function [SD,labels,bounds] = segmentCellMonolayer(filePath,varargin)
%SEGMENTCELLMONOLAYER segments the cell monolayer image given in filename
%   using parameters of gaussian blur amount, blob deletion size, 
%   erosion/dilation shape size, and reconstruction shape size
%
%   **REQUIRES MATLAB IMAGE PROCESSING TOOLBOX**
%
%INPUTS: 
%   filePath:   Filepath of folder containing TIF images to be segmented 
%       specified as a character array. All images contain in the specified 
%       folder will be analyzed unless filerange is specified using startIm 
%       and endIm. Images are analyzed and identified in ASCII dictionary 
%       order (be careful when numbering without leading zeros). 
%       (Use single quotes around this variable as it is a char argument).
%
%   startIm:    Specifies the number of the first image to be analyzed.
%       Files are numbered according to ASCII dictionary order.
%
%   endIm:      Specifies the number of the last image to be analyzed.
%       Files are numbered according to ASCII dictionary order.
%
%   frameSkip:  Specifies the number of images skipped inbetween image 
%       analyses. (i.e. every Nth image is analyzed)
%
%   blurAmt:    Adjusts the amount of gaussian blurring applied. Can 
%       improve segmentation by reducing noise and edge roughness
%       (default = 5)
%
%   blobSize:   Blob threshold size, in pixels, below which blobs are 
%       deleted. Can improve segmentation by removing unnecessary small
%       areas, preventing them from being detected as unique cells.
%       (default = 50)
%
%   gapAmt:     Controls the amount of erosion/dilation performed during 
%       the gap closing step of segmentation. Can improve segmentation by
%       reconnecting cell edges that were separated during initial
%       thresholding. (default = 2)
%   
%   flatAmt:    Controls the amount of erosion/dilation performed during 
%       the opening and closing-by-reconstruction steps which highlight and
%       flatten foreground objects. Can improve segmentation by altering
%       the level of flattening applied. (default = 4)
%
%   finalDil:   Controls the amount of dilation applied to the flattened
%       foreground objects, effectively smooths edges but does not greatly
%       affect segmentation. (default = 5)
%
%   boundColor: Controls the color of the boundaries in the segmented
%       image. Accepts the following example specifications: 'red', 'r', or
%       [1 0 0]. (default = 'green')
%
%   doPlot:     Controls whether the segmented images will also be 
%       displayed in Matlab figure windows in addition to being saved to
%       .tif files within a subfolder in the directory selected using
%       filePath. (useful for testing to quickly visualize segmentation
%       results for small test batches, use carefully with large image sets
%       as the large number of figures can cause hardware issues.)
%       (default = 0, no figure window creation)
%
%   imWriting:  Controls whether images are written to TIF files which are
%       placed in the specified directory and named according to the 
%       analyzed image appended with either '_labels' or '_bounds'.
%       (default = 0, no image writing)
%
%   stackMode:  Controls whether images are written to a single multipage
%       TIF file or multiple TIF files (one per analyzed frame). 
%       (default = 0, writes images to multiple TIF files)
%
%   outRes:     Controls the output resolution of the final label and
%   boundary segmented images. (default = 600)
%
%
%OUTPUTS:
%   SD:         Segmentation data structure that contains the labels and
%       bounds variables, the structure is also populated further in
%       downstream functions and serves to transfer data between functions.
%
%   labels:     Cell array of segmented cells overlaid with unique colors
%       identifying each individual cell. The labels output is designed for
%       use with the next function in the monolayer analysis series
%       (calculateMonolayerStatistics.m) This output is also written to 
%       SD.mat and can be replaced with '~' in the function call if not
%       needed in the workspace.
%
%   bounds:     Cell array of cell images overlaid with segemntation
%       boundaries. This output is also written to SD.mat and can be 
%       replaced with '~' in the function call if not needed in the 
%       workspace.
%
%EXAMPLE:
%   [SD,labels,bounds] = segmentCellMonolayer('/storage/disk1/sehaarma/GP065 SET',1,7,2,~,100,~,5,4,'red',0,1,0,~)
%       Analyzing images 1 to 7 skipping every other image (4 total
%       analyses), using default blur amount, 100 pixel blob thresholding,
%       default gap closing, increased flattening, decreased final
%       smoothing, red boundary color, plotting off, image writing on,
%       image stacking off, and default output resolution.

%% ** SET DEFAULT VALUES FOR OPTIONAL INPUTS ** 
imagefiles = dir(fullfile(filePath, '*.tif'));
nImages = numel(imagefiles);
defaultStartIm = 1;
defaultEndIm = nImages;
defaultFrameSkip = 1;
defaultBlurAmt = 5;
defaultBlobSize = 50;
defaultGapAmt = 2;
defaultFlatAmt = 4;
defaultFinalDil = 5;
defaultBoundColor = 'green';
defaultDoPlot = 0;
defaultImWriting = 0;
defaultStackMode = 0;
defaultOutRes = 600;

%% ** PARSE INPUTS ** 
disp('Running segmentCellMonolayer')
fprintf('Parsing Inputs...');
p = inputParser;
validBoundColor = @(x) isstring(x) || isnumeric(x) || ischar(x);
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
validScalarPosNumNonMax = @(x) isnumeric(x) && isscalar(x) && (x > 0) && (x <= nImages);
validBinary = @(x) (x == 1) || (x == 0);
addRequired(p,'filePath',@(x)exist(x,'dir'));
addOptional(p,'startIm',defaultStartIm,validScalarPosNum);
addOptional(p,'endIm',defaultEndIm,validScalarPosNumNonMax);
addOptional(p,'frameSkip',defaultFrameSkip,validScalarPosNumNonMax);
addOptional(p,'blurAmt',defaultBlurAmt,validScalarPosNum);
addOptional(p,'blobSize',defaultBlobSize,validScalarPosNum);
addOptional(p,'gapAmt',defaultGapAmt,validScalarPosNum);
addOptional(p,'flatAmt',defaultFlatAmt,validScalarPosNum);
addOptional(p,'finalDil',defaultFinalDil,validScalarPosNum);
addOptional(p,'boundColor',defaultBoundColor,validBoundColor);
addOptional(p,'doPlot',defaultDoPlot,validBinary);
addOptional(p,'imWriting',defaultImWriting,validBinary);
addOptional(p,'stackMode',defaultStackMode,validBinary);
addOptional(p,'outRes',defaultOutRes,validScalarPosNum);
parse(p,filePath,varargin{:});

startIm = p.Results.startIm;
endIm = p.Results.endIm;
frameSkip = p.Results.frameSkip;
blurAmt = p.Results.blurAmt;
blobSize = p.Results.blobSize;
gapAmt = p.Results.gapAmt;
flatAmt = p.Results.flatAmt;
finalDil = p.Results.finalDil;
boundColor = p.Results.boundColor;
doPlot = p.Results.doPlot;
imWriting = p.Results.imWriting;
stackMode = p.Results.stackMode;
outRes = p.Results.outRes;
disp(' done!');

%% ** SEGMENT IMAGES **
disp('Beginning image segmentation...');
%Initialize label and boundary cell arrays
frameList = startIm:frameSkip:endIm; %generate numerical list of frames
totalFrames = length(frameList); %extract total # of frames to be analyzed
labels{totalFrames} = [];
bounds{totalFrames} = [];
%Initialize command window writing to display the image number being analyzed
reverseStr = ''; fmt_length = num2str(numel(num2str(endIm))); fmt_spec = ['%' fmt_length '.0f'];
for ii = 1:totalFrames
    if ismember(ii,1)
        blobSize = 80;
        flatAmt = 7;
        finalDil = 6;
        blurAmt = 5;
    elseif ismember(ii,[5,9,10,11,12,13])
        blobSize = 200;
        flatAmt = 9;
        finalDil = 9;
        blurAmt = 7;
    else
        blobSize = 80;
        blurAmt = 5;
        gapAmt = 5;
        flatAmt = 8;
        finalDil = 6;
    end
    %Write current image number to command window with readable formatting
    msg = sprintf([' frame number: ' fmt_spec '/' num2str(length(frameList))], ii);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));

    currFileName = imagefiles(frameList(ii)).name; %get image filename from framelist
    currIm = imread(currFileName); %load image
    currIm = rgb2gray(currIm);

    %Perform initial thresholding to generate background markers
    [histCounts,~] = imhist(currIm,256); %get histogram of unaltered pixel intensities
    threshold = otsuthresh(histCounts); %convert intensity histogram to threshold

    adjIm = imadjust(currIm); %optimize image contrast
    blurIm = imgaussfilt(adjIm,blurAmt); %apply gaussian blurring

    bwIm = imbinarize(blurIm,threshold); %thresholding

    %Perform blob removal via image opening
    bwIm = bwareaopen(bwIm,blobSize); %image opening

    %Perform gap closing via erosion/dilation to improve background markers
    dilShape = strel('disk',gapAmt);
    eroShape = strel('disk',gapAmt-1);
    bwIm = imdilate(bwIm,dilShape); %dilation
    bwIm = imerode(bwIm,eroShape); %erosion
    bwIm = imdilate(bwIm,dilShape); %dilation
    bwIm = imerode(bwIm,eroShape); %erosion

    %Perform flattening of image local minima to generate markers
    shape = strel('disk',flatAmt);
    eroIm = imerode(blurIm,shape);
    recIm = imreconstruct(eroIm,blurIm);
    dilIm = imdilate(recIm,shape);
    flatIm = imreconstruct(imcomplement(dilIm),imcomplement(recIm));
    flatIm = imcomplement(flatIm);

    %Generate foreground markers and constrain them to desired locations
    fgmIm = imregionalmin(flatIm,4);
    fgmIm = imdilate(fgmIm,strel('disk',finalDil));
    prewsIm = imimposemin(double(bwIm),double(fgmIm));

    %Perform watershed segmentation
    wslabels = watershed(prewsIm);

    %Obtain label matrix
    labels{ii} = wslabels;

    %Obtain boundary matrix
    boundmask = boundarymask(wslabels);
    bounds{ii} = boundmask;

    if doPlot
        figure, imshow(labeloverlay(adjIm,wslabels))
        figure, imshow(imoverlay(adjIm,boundmask,boundColor))
    end
    
    if ~isfolder(strcat(filePath,filesep,'segmentedImages'))
        mkdir(strcat(filePath,filesep,'segmentedImages'));
    end

    if ~stackMode && imWriting
        imwrite(labeloverlay(adjIm,wslabels),strcat(filePath,filesep,'segmentedImages',filesep,currFileName(1:end-4),'_labels.tif'),'Resolution',outRes)
        imwrite(imoverlay(adjIm,boundmask,boundColor),strcat(filePath,filesep,'segmentedImages',filesep,currFileName(1:end-4),'_bounds.tif'),'Resolution',outRes)
        %imwrite(labels{ii},strcat(filePath,filesep,'segmentedImages',filesep,currFileName(1:end-4),'_labels.tif'),'Resolution',outRes)
        %imwrite(bounds{ii},strcat(filePath,filesep,'segmentedImages',filesep,currFileName(1:end-4),'_bounds.tif'),'Resolution',outRes)
    end
end
disp(' done!');

%% ** SAVE RESULTS **
fprintf('Saving results...');
SD.labels = labels; SD.bounds = bounds;
save(strcat(filePath,filesep,'SD.mat'),'-struct','SD');
disp(' done!');

if stackMode && imWriting
    imwrite(labels{1},strcat(imagefiles(frameList(1)).name(1:end-4),'_labels.tiff'),'Resolution',outRes)
    imwrite(bounds{1},strcat(imagefiles(frameList(1)).name(1:end-4),'_bounds.tiff'),'Resolution',outRes)
    for jj = 2:totalFrames
        imwrite(labels{jj},strcat(imagefiles(frameList(jj)).name(1:end-4),'_labels.tiff'),'Resolution',outRes,'WriteMode','append')
        imwrite(bounds{jj},strcat(imagefiles(frameList(jj)).name(1:end-4),'_bounds.tiff'),'Resolution',outRes,'WriteMode','append')
    end
end
end