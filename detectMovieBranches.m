function  movieData = detectMovieBranches(movieData,paramsIn)













if ~isa(movieData,'MovieData3D')
    error('The first input argument must be a valid MovieData object for a 3D movie!')
end

%Make sure the movie has masks
iSegProc = cellfun(@(x)(isa(x,'SegmentationProcess3D')),movieData.processes_);

if isempty(iSegProc)
    error('The movie must be segmented before the branches can be found!')
end

% Check if the branch-detection has been run before
% iProc = cellfun(@(x)(isa(x,'BranchDetectionProcess')),movieData.processes_);
% 
% if isempty(iProc)
%     iProc = numel(movieData.processes_)+1;
%     movieData.addProcess(BranchDetectionProcess(movieData));
% end
% 
% 
% if nargin < 2
%     paramsIn = [];
% end
% 
% p = parseProcessParams(movieData.processes_{iSegProc},paramsIn);

maskNames = movieData.processes_{iSegProc}.getMaskFileNames(1);
maskDir = movieData.processes_{iSegProc}.maskPaths_{1};
nFrames = movieData.nFrames_;

MakeQTMovie('start',[movieData.outputDirectory_ filesep 'branch movie.mov']);

for j = 1:nFrames
    
    mask = tif3Dread([maskDir filesep maskNames{1}{j}]);
    
    branches(j) = getBranchesFromMask(mask);    
    view([-38.5,42]);
    MakeQTMovie('addfigure');
    
end

MakeQTMovie('finish');

save([movieData.outputDirectory_ filesep 'branch detection.mat'],'branches');
