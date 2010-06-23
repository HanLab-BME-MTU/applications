function make3DMaskMovie(movieData,iChan)
%MAKE3DMASKMOVIE makes a movie overlaying the 3D masks on the fluorescence images
% 
% make3DMaskMovie(movieData,iChan)
% 
% This function makes a movie overlaying the segmentation on
% the fluorescence for the input movie and channel.
% 
% Input:
% 
%   movieData - A MovieData3D object describing the movie.
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

if ~isa(movieData,'MovieData3D')
    error('The first input must be a valid MovieData3D object!')
end

if nargin < 2 || isempty(iChan)
    iChan = 1;
end

%Check for existing seg processes
iSegProc = cellfun(@(x)(isa(x,'SegmentationProcess3D')),movieData.processes_);

if isempty(iSegProc)
    error('The movie has not been segmented! Please segment movie first!')
end

if ~movieData.processes_{iSegProc}.checkChannelOutput(iChan);
    error('The selected channel does not have valid masks! Check specified channel!')
end

%% ----- Init ----- %%


%Get the mask and image file names
maskFiles = movieData.processes_{iSegProc}.getMaskFileNames(iChan);
maskDir = movieData.processes_{iSegProc}.maskPaths_{iChan};
imageFiles = movieData.getImageFileNames(iChan);
imageDir = movieData.channelPath_{iChan};


%% ---- Make the movie ----- %%

nImages = movieData.nFrames_;


for iImage = 1:nImages;
   
    clf    
    
    %Load image & mask for this timepoint
    currIm = double(stackRead([imageDir filesep imageFiles{1}{iImage}]));    
    currMask = tif3Dread([maskDir filesep maskFiles{1}{iImage}]);
    
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
    maskCaps = isocaps(currMask,.5);
    patch(maskCaps,'EdgeColor','none','FaceColor','b','FaceAlpha',.5)
    light %Add a light    
    
    title(['Frame ' num2str(iImage)])
    
    if iImage == 1
        MakeQTMovie('start',[movieData.outputDirectory_ filesep mvName '.mov'])
        MakeQTMovie('quality',.9)                    
        
    end
    
    MakeQTMovie('addfigure')
    
end

MakeQTMovie('finish')


