function movieData = analyze3DMovieMaskGeometry(movieData,paramsIn)
%ANALYZE3DMASKGEOMETRY calculates various properties of the 3D masks for the input movie
% 
% movieData3D = analyze3DMovieMaskGeometry(movieData3D);
% movieData3D = analyze3DMovieMaskGeometry(movieData3D,paramIn);
% 
% This function calls analyze3DMaskGeometry.m on every mask of the input
% movie and writes the output to files. See analyze3DMaskGeometry.m for
% details.
% 
% *NOTE:* If the MovieData has specified xy and z pixel sizes, the mask
% properties will take into account the voxel aspect ratio. If not
% available, the pixel sizes and z spacing will be assumed to be 1. This
% may cause innacurate geometric properties.
%
% Input:
%   
%   movieData3D - A MovieData3D objec describing the movie to analyze the
%   masks from. This move must have already been segmented with
%   segment3DMovie.m
% 
%   paramsIn - A structure containing optional parameters. The possible
%   paramter field names and values are:
% 
%       (FieldName->fieldValue)
%
%       ('ChannelIndex' -> positive integer) Integer index of the
%       channel to analyze masks from.
%       Optional. Default is channel 1.
%
%       ('PhysicalUnits'->true/false) If true, and the movie has pixel
%       sizes and z-spacing specified, then the geometric properties will
%       be specified in real-word physical units (nm, nm^3, 1/nm etc.) Note
%       that even if this is  set to false, if the movie has pixel size and
%       z-spacing specified, the properties will still take into account
%       the voxel asymmetry. This means that the coordinates will not be in
%       matrix coordinates.
%       Optional. Default is false.
%
%       ('OutputDirectory' -> character string) Optional. A character
%       string specifying the directory to save the analysis to.
%       Optional. Default is a sub-directory of the movie's outputDirectory
%       called "mask_geometry"
%
%       (BatchMode->true/false) If true, all graphical output is
%       suppressed, such as progress bars, figures etc..
%
%       ***NOTE:*** Additionally, all of the parameters used by
%       analyze3DMaskGeometry.m can be specified as fields, where the field
%       names are the same as the option names in that function.
% 
% Output:
% 
% movieData3D - The updated MovieData3D object, with the analysis logged in
% the processes_ array.
% 
% Additionally, the analysis results will be written to file in the
% location specified by the OutputDirectory parameter.
% 
% Hunter Elliott
% 3/2011
% 

%% ------------------------- Parameters -------------------------%%

fName = 'mask_geometry_frame_';%String for naming geometry files.



%% ---------------------- Input ------------------------------- %%

if nargin < 1 || ~isa(movieData,'MovieData3D')
    error('The first input must be a valid MovieData3D object!')
end

if nargin < 2 || isempty(paramsIn)
    paramsIn = [];
end

%Get the indices of any previous processes from this function
iProc = movieData.getProcessIndex('MaskGeometry3DProcess',1,0);

%If the process doesn't exist, create it with default settings.
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(MaskGeometry3DProcess(movieData));                                                                                                 
end

%Parse input, store in parameter structure
p = parseProcessParams(movieData.processes_{iProc},paramsIn);

if numel(p.ChannelIndex) > 1
    error('You may only specify one channel!')
end



%% ---------------------- Init ------------------------------%%


%Make sure the movie has been segmented and has masks
iSegProc = movieData.getProcessIndex('SegmentationProcess3D',1,~p.BatchMode);

if isempty(iSegProc) || ~movieData.processes_{iSegProc}.checkChannelOutput(p.ChannelIndex)
    error('The input MovieData does not have valid masks for the selected channel!');
end

%Get mask file names and directory
maskDir = movieData.processes_{iSegProc}.outFilePaths_{1,p.ChannelIndex};
maskNames = movieData.processes_{iSegProc}.getOutMaskFileNames(p.ChannelIndex);

nFrames = movieData.nFrames_;

if ~p.BatchMode
    wtBar = waitbar(0,'Please wait, analyzing mask geometry...');
end        

%Format string for zero-padding file names
fString = ['%0' num2str(floor(log10(nFrames))+1) '.f'];

%Get the xy and z pixel sizes for scaling
if ~isempty(movieData.pixelSize_) && ~isempty(movieData.zSpacing_)
    pixXY = movieData.pixelSize_;
    pixZ = movieData.zSpacing_;
    hasSizes = true;    
else
    %warn the user, and assume unit pixel sizes.
    warning('Migration3D:MissingVoxelDimensions',...
        'Pixel XY size and Z spacing not specified in input movieData! Geometry analysis will assume symmetric voxels of size 1nm!');
    pixXY = 1;
    pixZ = 1;
    hasSizes = false;
end


outDir = p.OutputDirectory;
mkClrDir(outDir);

disp('Starting mask geometry analysis...')
disp(['Using masks from ' maskDir ', results will be saved to ' outDir ])


%% ---------------------- Analysis --------------------------%%


for iFrame = 1:nFrames
    
    %Load the current mask
    currMask = tif3Dread([maskDir filesep maskNames{1}{iFrame}]);
                
    
    if pixXY ~= pixZ
        %Scale the mask so that the voxel aspect ratio is taken into account
        currMask = make3DImageVoxelsSymmetric(currMask,pixXY,pixZ);
    end
    %Get the geometry properties
    maskProp = analyze3DMaskGeometry(currMask,'SmoothSigma',p.SmoothSigma,...
                        'IsoValue',p.IsoValue); %#ok<*NASGU>
    
    if hasSizes && p.PhysicalUnits
        %Scale these properties so that they are in nm. We have scaled the
        %matrix in the z direction to get symmetric voxels now so we scale
        %all dimensions by the xy pixel size.
        nObj = numel(maskProp);
        for i = 1:nObj
            maskProp(i).SmoothedSurface.vertices = maskProp(i).SmoothedSurface.vertices .* pixXY;
            maskProp(i).SurfaceNorms = maskProp(i).SurfaceNorms .* pixXY;
            maskProp(i).GaussianCurvature = maskProp(i).GaussianCurvature ./ (pixXY^2);
            maskProp(i).CurvaturePC1 = maskProp(i).CurvaturePC1 ./ pixXY;
            maskProp(i).CurvaturePC2 = maskProp(i).CurvaturePC2 ./ pixXY;
            maskProp(i).Centroid = maskProp(i).Centroid .* pixXY;
            maskProp(i).Volume = maskProp(i).Volume * pixXY^3;
            maskProp(i).CenterMostDist = maskProp(i).CenterMostDist * pixXY;
            maskProp(i).CenterMost = maskProp(i).CenterMost .* pixXY;
            maskProp(i).MeanCurvature = maskProp(i).MeanCurvature ./ pixXY;
        end
    end
    
    
    
    %Save them to file
    numStr = num2str(iFrame,fString);
    save([outDir filesep fName numStr '.mat'],'maskProp');    
  
    if ~p.BatchMode
        %Update the waitbar occasionally to minimize slowdown
        waitbar(iFrame/nFrames,wtBar)
    end
        
    
end

%% ---------- Finalization ---------- %%

%Store the input/output directories in the movieData
movieData.processes_{iProc}.setOutFilePath(p.ChannelIndex,outDir);
movieData.processes_{iProc}.setInImagePath(p.ChannelIndex,maskDir);

if ~p.BatchMode && ishandle(wtBar)
    close(wtBar)
end

movieData.processes_{iProc}.setDateTime;
movieData.saveMovieData;

disp('Finished analyzing mask geometry!')






