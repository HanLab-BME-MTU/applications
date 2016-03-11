%% program to rotate diaSLM data. Additionally can interpolate data in z to a finer grid for low NA mode
%
%
%Reto Fiolka, 07-09-2015
%
%Attempted to turn into a function, with automatic calling of image
%features.  Kevin Dean, 07-29-2015.
%
%inDirectory is the directory that the files exist in, terminated with a
%slash.
%
%interpolate = 1 if yes, 0 if no.
%
%Good start for xDeg=44, yDeg=2.
%
%If test = 1, then it will only process the first image.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = diagonolAffine3D(inDirectory,interpolate,zStep,xDeg,yDeg,test)

files = dir(inDirectory); % find the names of the files in the directory

savePath = [inDirectory filesep 'affineTransformed' filesep]; % set the save directory

if ~isdir(savePath) % make a directory to save images in
    mkdir([savePath])
end

imageIndex = 1; % find the names of the images

imageNames = [];

for i = 1:length(files) % check every file in the directory
    if strncmp(fliplr(files(i).name), 'fit.', 4);% if the image name ends in '.tif', save the name
        imageNames{imageIndex} = files(i).name;
        imageIndex = imageIndex + 1;
    end
end

warning('off', 'MATLAB:imagesci:tifftagsread:expectedAsciiDataFormat');


for imageNumber = 1:length(imageNames)
    
    
    % display progress
    if mod(imageNumber,2) == 0
        disp(['  Reading ' num2str(imageNumber) '/' num2str(length(imageNames))])
    end
    
    % find the image size
    imageInfo = imfinfo([inDirectory imageNames{imageNumber}]);
    imageSize = [imageInfo(1).Height; imageInfo(1).Width; length(imageInfo)];
    
    % initiate the image variable
    image3D = zeros(imageSize(1), imageSize(2), imageSize(3));
    
    % load each plane in Z
    parfor z = 1:imageSize(3)
        image3D(:,:,z) = im2double(imread([inDirectory imageNames{imageNumber}], 'Index', z, 'Info', imageInfo));
    end
    
    %% Interpolate Z-axis
    
    if interpolate == 1;
        
        zInterpolate=round(imageSize(3)/0.160*zStep);
        
        image3DInterpolated=zeros(imageSize(1),imageSize(2),zInterpolate);
        
        parfor i=1:imageSize(2)
            
            tempG=squeeze((image3D(:,i,:)));
            
            image3DInterpolated(:,i,:)=imresize(tempG,[imageSize(1),zInterpolate]);
            
        end
        
        image3D=image3DInterpolated;
        
        clear image3DInterpolated
        
    end
    
    %% rotate axis
    
    xRot = -xDeg*pi/180;         % main rotation of ~45 degrees
    
    yRot = -yDeg*pi/180;       % small rotation for sample mount tilt
    
    rotationMatrixX = [1 0 0 0; 0 cos(xRot) -sin(xRot) 0; 0 sin(xRot) cos(xRot) 0; 0 0 0 1];
    
    rotationMatrixY=  [cos(yRot) -sin(yRot) 0 0; sin(yRot) cos(yRot) 0 0; 0 0 1 0; 0 0 0 1];
    
    rotationMatrix=rotationMatrixX*rotationMatrixY;
    
    tform = affine3d(rotationMatrix);
    
    outputImage=imwarp(image3D,tform);
    
    outputImage=permute(outputImage,[3,2,1]);
    
    %%truncate empty space above and below cell
    midZ=round(size(outputImage,3)/2);
    
    outputImage=outputImage(:,:,midZ-50:midZ+50);
    
    %% Save the Affine Transformed Data
    
    savePath = [inDirectory filesep 'affineTransformed' filesep];
    
    stackCounter = 0;
    
    for i = 1:1:size(outputImage,3);
        
        stackCounter=stackCounter+1;
        
        tiffStructure=Tiff([savePath 'transformed_' imageNames{1}],'a');
        
        im2Save=im2uint16(outputImage(:,:,i));  %im2uint16 command for 16-bit image.
        
        tagstruct.ImageLength = size(outputImage,1);
        
        tagstruct.ImageWidth = size(outputImage,2);
        
        tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
        
        tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
        
        tagstruct.Software = 'MATLAB';
        
        tagstruct.BitsPerSample = 16;
        
        tiffStructure.setTag(tagstruct);
        
        tiffStructure.write(im2Save);
        
        tiffStructure.close();
        
    end
    
    %% Save the Affine Transformed MIP Data
    
    maxXY = squeeze(max(outputImage, [], 3));
    
    % set the directory where images will be saved (you can change this)
    savePath = [inDirectory filesep 'affineTransformedMIP' filesep];
    
    % make a directory to save images in
    if ~isdir(savePath)
        mkdir([savePath])
    end
    
    tiffStructure=Tiff([savePath 'transformedMIP_' imageNames{1}],'a');
    
    im2Save=im2uint16(maxXY(:,:));  %im2uint16 command for 16-bit image.
    
    tagstruct.ImageLength = size(outputImage,1);
    
    tagstruct.ImageWidth = size(outputImage,2);
    
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    
    tagstruct.Software = 'MATLAB';
    
    tagstruct.BitsPerSample = 16;
    
    tiffStructure.setTag(tagstruct);
    
    tiffStructure.write(im2Save);
    
    tiffStructure.close();
    
    if test == 1
        break
    end
    
end

% turn the warning back on
warning('on','MATLAB:imagesci:tifftagsread:expectedAsciiDataFormat');