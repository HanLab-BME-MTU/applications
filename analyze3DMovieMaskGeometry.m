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
% *NOTE:* The mask property structures which are written to disk do not
% take into account the voxel size or aspect ratio. That is, any x-y-z
% coordinates will be stored as if the xy and z pixel sizes were all equal
% to 1.
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
maskDir = movieData.processes_{iSegProc}.outMaskPaths_{p.ChannelIndex};
maskNames = movieData.processes_{iSegProc}.getOutMaskFileNames(p.ChannelIndex);

nFrames = movieData.nFrames_;

if ~p.BatchMode
    wtBar = waitbar(0,'Please wait, analyzing mask geometry...');
end        

%Format string for zero-padding file names
fString = ['%0' num2str(floor(log10(nFrames))+1) '.f'];

outDir = p.OutputDirectory;
mkClrDir(outDir);

disp('Starting mask geometry analysis...')
disp(['Using masks from ' maskDir ', results will be saved to ' outDir ])

%%------------------------- Parameters -------------------------%%

fName = 'mask_geometry_frame_';%String for naming geometry files.


%% ---------------------- Analysis --------------------------%%


for iFrame = 1:nFrames
    
    %Load the current mask
    currMask = tif3Dread([maskDir filesep maskNames{1}{iFrame}]);
                
    %Get the geometry properties
    maskProp = analyze3DMaskGeometry(currMask,'SmoothSigma',p.SmoothSigma,...
                        'IsoValue',p.IsoValue); %#ok<*NASGU>
    
    %Save them to file
    numStr = num2str(iFrame,fString);
    save([outDir filesep fName numStr '.mat'],'maskProp');    
  
    if ~p.BatchMode && mod(iFrame,5)
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






