function deconvolveDirectory(inDirectory, savePath)

% deconvolveDirectory weiner deconvolves the images in the specified directory

% INPUTS:
%
% inDirectory - The directory in which images are saved. Images are assumed
%               to be tiff files with each file a 3D image at one
%               timepoint.
%
% savePath - The directory where deconvolved images are saved.  Images are 
%            saved with the same name as the orginal and as 16-bit tiff 
%            images.
% 
% NOTE that the PSF must be saved in a .mat file at the path specified by 
% PSFpath.  The PSF is assumed to be a variable named PSF whose dimensions
% are in the order (x,y,z). maxOTF, the second largest OTF value in the fft 
% of the full size PSF, is also assumed to be saved in this mat file. 
%
% Written by Reto Fiolka and revised by Meghan Driscoll (August 2014).


%% Set Parameters
PSFpath = '/project/cellbiology/gdanuser/melanomaModel/Meghan/Deconvolved/PSF/rotAvgPSF.mat'; % (Ses the NOTE above.)
measureWeiner = 1; % if 1 estimate the weiner parameter (this should only be enabled for mostly uniformly bright cells)

% deconvolution parameters
apoHeight = 0.04; %0.06, 0.012 the height of the apodization mask as a percentage of the OTF maximum (excluding the peak at 0 frequency)
weiner = 0.01; % the weiner parameters (only used if measureWeiner is set to 0)

% (You probably don't want to change the following paramaters)
% segmentation parameters (only used if measureWeiner is enabled)
morphRadiusErode = 5; % shrinks the estimated boundary by this number of pixels before measuring the mean signal
morphRadiusDilate = 20; % expands the estimated boundary by this number of pixels before measuring the noise variance

% hacked parameters
imageMultiply = 200;  % multiply the image by this number before saving it


%% Find the names of the images in the directory 
% (This is not yet MovieData compatible.)
disp('Finding images to deconvolve')

% make a directory in which to save deconvolved images
if ~isdir(savePath)
    mkdir(savePath)
end
    
% find the names of the files in the directory
files = dir(inDirectory);

% find the names of the images
imageIndex = 1;
imageNames = [];
for i = 1:length(files) % check every file in the directory
    
    % if the image name ends in '.tif', save the name
    if strncmp(fliplr(files(i).name), 'fit.', 4)  
        imageNames{imageIndex} = files(i).name;
        imageIndex = imageIndex + 1;
    end  

end


%% Find the Image Sizes
% All images in the directory are assumed to have the same size as the first image, since this assumption makes the code much faster.
imageInfo = imfinfo([inDirectory imageNames{1}]);
imageSize = [imageInfo(1).Height; imageInfo(1).Width; length(imageInfo)];


%% Load the PSF
disp('Loading the PSF')

% load the point spread function (PSF)
load(PSFpath);

% check that the PSF path contains a varible named PSF
if isempty(who('PSF'))
    error(['PSF must be a Matlab variable saved in ' PSFpath])
end

% check that the PSF path contains a varible named maxOTF
if isempty(who('maxOTF'))
    error(['maxOTF must be a Matlab variable saved in ' PSFpath])
end

% The saved PSF must be the same size or larger than the images
sizePSF = size(PSF);
if (sizePSF(1)<imageSize(1)) || (sizePSF(2)<imageSize(2)) || (sizePSF(3)<imageSize(3))
    error('The PSF must the same size or larger than the images.')
end


%% Calculating the OTF
disp('Calculating the OTF')

% find the box of size imageSize about the PSF origin
originPSF = ceil((sizePSF+ones(1,3))/2); % (for even sizes the origin occurs above the center, i.e. the origin of an image with size 4x4 occurs at (3,3) )
PSF = PSF-min(PSF(:));
PSF = PSF./max(PSF(:));
smallPSF = PSF((originPSF(1)-ceil((imageSize(1)-1)/2)):(originPSF(1)+floor((imageSize(1)-1)/2)), ...
    (originPSF(2)-ceil((imageSize(2)-1)/2)):(originPSF(2)+floor((imageSize(2)-1)/2)), ...
    (originPSF(3)-ceil((imageSize(3)-1)/2)):(originPSF(3)+floor((imageSize(3)-1)/2))); clear PSF;

% find the OTF
OTF = fftshift(fftn(smallPSF)); clear smallPSF ;
OTF = abs(OTF);


%% Create an apodizationan filter
% create a triangular apodization filter that decays to zero at the edge of a window defined by a height threshold of the OTF
disp('Making the apodization mask')

% smooth the OTF
smoothedOTF = filterGauss3D(OTF,3); % (the Gaussian sigma is hardcoded, but the goal of this smoothing is simply to dampen the artifacts along the axes)

% make a binary apodization mask
apodizeBinary = smoothedOTF>(maxOTF*apoHeight); clear smoothedOTF;

% find the distance from the center to each pixel
centerDist = zeros(imageSize(1), imageSize(2), imageSize(3));
centerDist(ceil((imageSize(1)+1)/2), ceil((imageSize(2)+1)/2), ceil((imageSize(3)+1)/2)) = 1;
centerDist = bwdist(centerDist);

% find the distance from the binarized apodization mask edge to each pixel
edgeDist = bwdist(~apodizeBinary);

% make the apodization mask
apodizeMask = ones(imageSize(1), imageSize(2), imageSize(3)) - centerDist./(centerDist+edgeDist); clear centerDist edgeDist

% remove values outside of the binaryMask
apodizeMask(apodizeBinary == 0) = 0; clear apodizeBinary
 
% % Debug figures
% figure
% imagesc(apodizeMask(:,:,1+ceil((imageSize(3)+1)/2)))
% axis equal
% 
% figure
% imagesc(log(OTF(:,:,1+ceil((imageSize(3)+1)/2))))
% axis equal
% 
% figure
% imagesc(squeeze(apodizeMask(:,1+ceil((imageSize(2)+1)/2),:)))
% axis equal
% 
% figure
% imagesc(squeeze(log(OTF(:,1+ceil((imageSize(2)+1)/2),:))))
% axis equal


%% Weiner deconvolve all of the images
disp('Deconvolving')

% create 3D spheres as structuring element for later dilations and erosions
[x,y,z] = ndgrid(-morphRadiusErode:morphRadiusErode,-morphRadiusErode:morphRadiusErode,-morphRadiusErode:morphRadiusErode);
sphereErodeSE = ((x.*x+y.*y+z.*z)./morphRadiusErode^2)<1; clear x y z;

[x,y,z] = ndgrid(-morphRadiusDilate:morphRadiusDilate,-morphRadiusDilate:morphRadiusDilate,-morphRadiusDilate:morphRadiusDilate);
sphereDilateSE = ((x.*x+y.*y+z.*z)./morphRadiusDilate^2)<1; clear x y z;

% iterate through the images
weinerEstimateList = []; % assemble a list of the estimated Weiner parameters
for n = 1:length(imageNames)
    
    % display progress
    disp(['   ' imageNames{n}]);
    
    
    % try to find information about the image
    try
        imageInfo = imfinfo([inDirectory imageNames{n}]);
    catch 
        disp([imageNames{n} ' is not an image and will not be deconvolved.'])
        continue
    end
    
    % check that the image size is as expected
    thisImageSize = [imageInfo(1).Height; imageInfo(1).Width; length(imageInfo)];
    if max(imageSize ~= thisImageSize)
        disp([imageNames{n} ' has an unexpected size and will not be deconvolved.'])
        continue
    end
        
    % load the image into a single large matrix
    image3D = zeros(imageSize(1), imageSize(2), imageSize(3));
    for z = 1:imageSize(3) 
        image3D(:,:,z) = im2double(imread([inDirectory imageNames{n}], 'Index', z, 'Info', imageInfo));
    end
    
    % try to estimate the weiner parameter
    if measureWeiner == 1
        
        % blur the image
        imageBlured = filterGauss3D(image3D,1); % blur the image
        imageBlured = imageBlured-min(imageBlured(:));
        imageBlured = imageBlured./max(imageBlured(:));
        
        % threshold the image
        imageThresh = imageBlured>graythresh(imageBlured(:)); clear imageBlured; % Otsu threshhold the image;
        
        % remove small objects in the thresholded image
        imageThresh = bwareaopen(imageThresh, 100);
        
        % measure the mean signal intensity
        imageErodeMask = imerode(imageThresh, sphereErodeSE); % erode the images
        imageSignal = imageErodeMask.*image3D; % find the mean signal intensity
        meanSignal = mean(imageSignal(imageErodeMask)); clear imageSignal imageErodeMask;
        
        % measure the variance of the background
        imageDilateMask = imdilate(imageThresh, sphereDilateSE); clear imageThresh;
        imageBackground = (~imageDilateMask).*image3D; % find the std of the background pixels
        stdBackground = std(imageBackground(~imageDilateMask));  clear imageBackground imageDilateMask;
        
        % estimate the weiner parameter
        weinerEstimate = stdBackground/meanSignal;
        weinerEstimateList = [weinerEstimateList, weinerEstimate];
    else
        
        % set the weiner parameter if it's not beinge estimated from the image
        weinerEstimate = weiner; 
    end
    
    % Fourier transform the image
    image3D = fftshift(fftn(image3D));

    % perform the Weiner deconvolution
    image3D = image3D.*OTF./((OTF.*OTF)+weinerEstimate);
    
    % perform apodization
    image3D = image3D.*apodizeMask;

    % inverse Fourier transform back to the image domain
    image3D = ifftn(ifftshift(image3D));

    % take the absolute value of the image
    image3D = abs(image3D);

    % remove image values less than 0
    image3D = image3D.*(image3D>0);
    
    % this is a hack
    image3D = imageMultiply*image3D;
    
    % convert the image to a 16-bit integer
    image3D = uint16((2^16-1)*image3D);
    
    % display warnings if the image dynamic range is incorrect
    if max(image3D(:)) < 100
        disp(['The maximum image value, ' num2str(max(image3D(:))) ', is low.'])
    elseif max(image3D(:)) == 2^16-1
        disp('The image brightness is saturated.')
    end
    
    % save the image
    imwrite(squeeze(uint16(image3D(:,:,1))),fullfile(savePath,imageNames{n}),'Compression','none') % overwrite any existing image
    for z=2:imageSize(3)
        imwrite(squeeze(uint16(image3D(:,:,z))),fullfile(savePath,imageNames{n}),'Compression','none','WriteMode','append')
    end

end

% display characteristics of the Weiner parameter
if measureWeiner == 1
    disp(['  Mean Weiner parameter: ' num2str(mean(weinerEstimateList))])
    disp(['  Std Weiner parameter: ' num2str(std(weinerEstimateList))])
    disp(['  Range Weiner parameter: ' num2str(max(weinerEstimateList)-min(weinerEstimateList))])
end

% save the Weiner estimate lists
save(fullfile(savePath, 'weiner.mat'), 'weinerEstimateList');