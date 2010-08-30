function  movieData = detectMovieBranches(movieData,paramsIn)
%DETECTMOVIEBRANCHES detects branches in each frame of the input 3D movie using getBranchesFromMask., 
% 
% movieData = detectMovieBranches(movieData,paramsIn)
% 
% movieData = detectMovieBranches(movieData)
% 
% This function calls getBranchesFromMask on every frame of the input movie
% and writes the resulting branches to the movie's output directory.
% 
% Input:
% 
%   movieData - A MovieData3D object describing the movie to be processed.
%
%   paramsIn - Structure with inputs for optional parameters. The
%   parameters should be stored as fields in the structure, with the field
%   names and possible values as described below
% 
%   Possible Parameter Structure Field Names:
%       ('FieldName' -> possible values)
%
%       ('OutputDirectory' -> character string)
%       Optional. A character string specifying the directory to save the
%       detected branches to.
%       If not input, the files will be saved to the same directory as the
%       movieData, in a sub-directory called "branch_detection"
%
%       ('ChannelIndex' -> Positive integer scalar) Optional. The
%       integer index of the channel to use masks from for branch
%       detection. If not input, the first channel will be used.
%
%       Pruning Parameters - Defaults and descriptions are in
%       getBranchesFromMask.m
%
%               ('MaxRadius->Positive Integer) 
% 
%               ('SmoothSigma->Positive Scalar) 
% 
%               ('IsoVal'->Positive Integer) 
%
%       ('BatchMode' -> True/False)
%       If true, graphical output and user interaction is
%       supressed (i.e. progress bars, dialog and question boxes etc.)
%
%
% Output:
%
%   movieData - the updated MovieData object with the parameters, paths
%   etc. stored in it, in the field movieData.processes_.
%
%   The branches are written to the directory specified by the parameter
%   OuptuDirectory, as a .mat file.
%
%
% Hunter Elliott, 
% 6/2010
%
%% ----- Parameters ----- %%

fName = 'branches_frame_'; %String for naming branch files for each frame


%% ---- Input ---- %%

if ~isa(movieData,'MovieData3D')
    error('The first input argument must be a valid MovieData object for a 3D movie!')
end

%Check if the branch-detection has been run before
iProc = movieData.getProcessIndex('BranchDetectionProcess',1,0);
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(BranchDetectionProcess(movieData));
end

if nargin < 2
    paramsIn = [];
end

%Resolve input and default parameters
p = parseProcessParams(movieData.processes_{iProc},paramsIn);

%Make sure the movie has masks
iSegProc = movieData.getProcessIndex('SegmentationProcess3D',1,~p.BatchMode);

if isempty(iSegProc)
    error('The movie must be segmented before the branches can be found!')
end


%Make sure the selected channel has valid masks
if ~movieData.processes_{iSegProc}.checkChannelOutput(p.ChannelIndex)
    error('The seclected channel does not have valid masks! Please check segmentation!');
end


%% ------ Init ------ %%

maskNames = movieData.processes_{iSegProc}.getOutMaskFileNames(p.ChannelIndex);
maskNames = maskNames{1}; %For parallelization - "slices" maskNames
maskDir = movieData.processes_{iSegProc}.outMaskPaths_{p.ChannelIndex};
nFrames = movieData.nFrames_;

disp('Starting branch detection...')
if ~p.BatchMode
    wtBar = waitbar(0,'Please wait, detecting branches...');
end

branches = cell(1,nFrames);


parfor j = 1:nFrames
    
    mask = tif3Dread([maskDir filesep maskNames{j}]);
    
    branches{j} = getBranchesFromMask(mask,p.MaxRadius,p.SmoothSigma,p.IsoValue);             %#ok<PFOUS>
   
%     if ~p.BatchMode
%        
%         waitbar(j/nFrames,wtBar)
%         
%     end
    
    disp(['Finished frame ' num2str(j)]);
    
end


%% ----- Finish ----- %%

if ~p.BatchMode && ishandle(wtBar)
    close(wtBar);
end


save([movieData.outputDirectory_ filesep 'detected_branches.mat'],'branches');


movieData.processes_{iProc}.setDateTime;
movieData.saveMovieData; %Save the new movieData to disk


disp('Finished branch detection!')


