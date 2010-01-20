function make3DMaskMovie(movieData,iChan)

% 
% make3DMaskMovie(movieData,iChan)
% 
% This function makes a movie overlaying the segmentation on
% the fluorescence for the input movie and channel.
% 
% Input:
% 
%   movieData - A movieData structure describing a 3D movie as created with
%   setup3DMovieData.m
% 
%   iChan - The index of the channel to make the movie from.
%   Optional. Default is 1.
% 
% 
% Output:
%    
%   The resulting movie will be saved as a .mov file in the movie's
%   analysis directory.  
% 
% 
% Hunter Elliott
% 1/2010
%

%% ----- Parameters ----- %%

mvName = 'mask_overlay_movie'; % The file name to save the movie as

%% ------- Input ------ %%

%Init / validate the movieData
movieData = setup3DMovieData(movieData);

if nargin < 2 || isempty(iChan)
    iChan = 1;
end

%Make sure it is a 3D movie
if ~isMovie3D(movieData)
    error('This function can only be used with 3D movies! Please create the movieData with setup3dMovieData!')
end

%Make sure the masks have been created
if ~checkMovieMasks(movieData,1);
    error('Problem with masks in input movie! Check mask files and movieData!')
end


%% ----- Init ----- %%


%Get the mask and image file names
maskFiles = getMovieMaskFileNames(movieData,iChan);
imageFiles = getMovieImageFileNames(movieData,iChan);


%% ---- Make the movie ----- %%

nImages = movieData.nImages(iChan);


for iImage = 1:nImages;
   
    clf    
    
    %Load image & mask for this timepoint
    currIm = double(stackRead(imageFiles{1}{iImage}));    
    currMask = tif3Dread(maskFiles{1}{iImage});
    
    %Normalize the image
%     currIm = currIm - min(currIm(:));
%     currIm = currIm * round((2^16 ./ max(currIm(:))));
%     currIm = uint16(currIm);
%     
    
    
    %Display the image
    subplot(1,2,1)
    viewStack(currIm)            
    
    %set the alpha axis 
    %This formula is completely arbitrary, and just happens to work
    %okay with my data....!!
    alim([median(currIm(:)) max(currIm(:))])        
    
    subplot(1,2,2)
    viewStack(currIm)               
    alim([median(currIm(:)) max(currIm(:))])
    
    
    %Overlay the mask
    maskSurf = isosurface(currMask,.5); %Get the surface of the mask
    patch(maskSurf,'EdgeColor','none','FaceColor','r','FaceAlpha',.5)
    light %Add a light    
    
    title(['Frame ' num2str(iImage)])
    
    if iImage == 1
        MakeQTMovie('start',[movieData.analysisDirectory filesep mvName '.mov'])
        MakeQTMovie('quality',.9)                    
        
    end
    
    MakeQTMovie('addfigure')
    
end

MakeQTMovie('finish')


