function imarisApp = viewMovieBranchesImaris(movieData3D,iChannel)
%VIEWMOVIEBRANCHESIMARIS overlays the branch detection on images using imaris 
% 
% viewMovieBranchesImaris(movieData);
% imarisApp = viewMovieBranchesImaris(movieData);
% 
% This function displays the results of branch detection and
% skeletonization overlain on the images from the input movie for visual
% inspection/validation.
% 
% Input:
% 
%   movieData3D - The MovieData3D object describing the movie to view.
% 
%   iChannel - The index of the channel to view in imaris, and the channel
%   which branch detection was run on. This corresponds to the channel
%   object's location within the channel array in the MovieData. This
%   channel must have valid masks. Optional. Default is channel 1.
% 
% Output:
%   
%   imarisApp - The handle to the imaris application. If deleted, imaris
%   will close.
%
%   Additionally, Imaris will be opened and then the images and detected
%   branches are loaded and displayed. The images will be displayed in red,
%   while the detected branch tips will be yellow spheres. The masks will
%   also be shown as white.
%   
%Hunter Elliott
%10/2010
%
%% ------- Parameters ------ %%

spotRad = 3; %Radius of spots showing detected branch tips in pixels;


%% -------- Input -------- %%

if nargin < 1 || isempty(movieData3D) || ~isa(movieData3D,'MovieData3D')
    error('The first input must be a valid MovieData3D object!')
end

if nargin < 2 || isempty(iChannel)
    iChannel = 1;
end

%Make sure the movie has valid branch detection
iProc = movieData3D.getProcessIndex('BranchDetectionProcess');
if isempty(iProc)
    error('The input MovieData must have a valid BranchDetectionProcess. Please run branch detection!')
end

%% ---------- Init ---------- %%

%Display the images in imaris first
imarisApp = viewMovieMasksImaris(movieData3D,iChannel);

%Load the branch data.
branches = load('C:\Users\HE19\Documents\Local_Data\3dmigration\bob_caax2_set2c_pos1\detected_branches_xferredover.mat'); %TEMPORARY, obviously! HLE
branches = branches.branches;

nImages = movieData3D.nFrames_;

%Create spots structure for detected branches
imarisSpots = imarisApp.mFactory.CreateSpots;
imarisSpots.mName = 'Detected Branch Tips';
imarisSpots.SetColor(1,1,0,0);

%Get pixel sizes, in nm
pixXY = movieData3D.pixelSize_;
pixZ = movieData3D.zSpacing_;



%% ----- Display Overlays ----- %%


%Determine number of branches in each frame and total
nPerFrame = cellfun(@(x)(size(x.tipLocation,1)),branches);
nTotal = sum(nPerFrame);

%Initialize spot arrays
branchTipXYZ = zeros(nTotal,3);
tipRadii = ones(nTotal,1)*pixXY*spotRad;
timeIndices = zeros(nTotal,1);             
currInd = 1;

%Go through each frame and add the branch detection overlays
for iImage = 1:nImages
    
    %Create time indices for these branches
    timeIndices(currInd:currInd+nPerFrame(iImage)-1) = iImage-1;
    %Store XY positions of branch tips
    branchTipXYZ(currInd:currInd+nPerFrame(iImage)-1,1:2) = ...
                branches{iImage}.tipLocation(:,1:2) * pixXY;
    branchTipXYZ(currInd:currInd+nPerFrame(iImage)-1,3) = ...
                branches{iImage}.tipLocation(:,3) * pixZ;            
    
    currInd = currInd+nPerFrame(iImage);
end

imarisSpots.Set(branchTipXYZ,timeIndices,tipRadii);
imarisApp.mSurpassScene.AddChild(imarisSpots);



