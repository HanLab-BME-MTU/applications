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
iSegProc = movieData.getProcessIndex('SegmentationProcess3D',1,1);

if isempty(iSegProc)
    error('The movie has not been segmented! Please segment movie first!')
end

if ~movieData.processes_{iSegProc}.checkChannelOutput(iChan);
     error('The selected channel does not have valid masks! Check specified channel!')
end

%% ----- Init ----- %%


%Get the mask and image file names
maskFiles = movieData.processes_{iSegProc}.getOutMaskFileNames(iChan);
maskDir = movieData.processes_{iSegProc}.outMaskPaths_{iChan};
imageFiles = movieData.getImageFileNames(iChan);
imageDir = movieData.getChannelPaths(iChan);


%% ---- Make the movie ----- %%

nImages = movieData.nFrames_;

fig = fsFigure(.75);

for iImage = 1:nImages;
   
    clf    
    
    %Load image & mask for this timepoint
    currIm = double(stackRead([imageDir{1} filesep imageFiles{1}{iImage}]));    
    currMask = tif3Dread([maskDir filesep maskFiles{1}{iImage}]);
    
    
    %Sub-sample so I can graduate
%     currIm = currIm(:,[1:2:end],:);
%     currIm = currIm([1:2:end],:,:);
%     currMask = currMask(:,[1:2:end],:);
%     currMask = currMask([1:2:end],:,:);
    
    %Display the image    
    viewStack(currIm)            
    view(-10,60)
    %Leave axis on four bounding box but turn ticks off
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    set(gca,'ZTick',[]);

    %set the alpha axis to allow some display saturation
    %This formula is completely arbitrary, and just happens to work
    %okay with my data....!!
    %if iImage == 1;
    alim([median(currIm(:)) max(currIm(:))/1.5])        
    %end

    
    %Overlay the mask
    maskSurf = isosurface(currMask,.5); %Get the surface of the mask
    patch(maskSurf,'EdgeColor','r','FaceColor','none','EdgeAlpha',.025)
    maskCaps = isocaps(currMask,.5);
    patch(maskCaps,'EdgeColor','b','FaceColor','none','EdgeAlpha',.015)
    %light %Add a light    
    
    title(['Frame ' num2str(iImage)])
    
    if iImage == 1
        MakeQTMovie('start',[movieData.outputDirectory_ filesep mvName '.mov'])
        MakeQTMovie('quality',.95)                    
        
    end
    
    MakeQTMovie('addfigure')
    
end

MakeQTMovie('finish')
close(fig);

