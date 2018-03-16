function [] = checkSegmentationAreaML(ML)
%function [] = checkSegmentationAreaML(ML) go over each movie in ML and
%check from MaskRefinement the area of the mask area, sort out too small
%areas, and print the movieData information regarding MDs with small
%segmentations for future manual correction. (Going over each MD through
%Segmentation Viewer is pain.
% Sangyoon Han 2018 March

ML = MovieList.load(ML.getFullPath);
numMDs = numel(ML.movies_);
%% Segmentation information
segAreas = zeros(numMDs,1);
mdNames = cell(numMD,1);

for ii=1:numMD
    curMD = ML.movies_{ii};
    % Get the MaskRefinement process
    curMaskRefProc = curMD.getProcess(curMD.getProcessIndex('MaskRefinementProcess'));
    iValidChan = find(curMaskRefProc.checkChannelOutput,1);
    curMask = curMaskRefProc.loadChannelOutput(iValidChan,1);
    segAreas(ii) = sum(curMask(:));
    mdNames{ii} = curMD.getFullPath;
end

%% Get statistics and get the low area movies
thresArea = median(segAreas)-std(segAreas);
idxLowAreas = segAreas<thresArea;
disp('MDs with Low Segmentation Areas:')
disp(mdNames(idxLowAreas))
newML = MovieList(ML.movies_(idxLowAreas),fullfile(ML.movieListPath_, 'ML_lowSegAreas.mat'));
newML.save
disp(['The newML is saved as ' newML.getFullPath])
