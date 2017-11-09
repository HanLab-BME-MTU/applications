function [ output_args ] = GCAChoosePatchSizeManuallyMovieList(projList)
%GCAChoosePatchSizeManuallyMovieList
% Small movieList wrapper designed to test multiple local threshold patch sizes
% and select a good candidate.
% Currently inputs a projList (simple cell not an object- eventually make for movieList Object)
% Assumes the LocalThreshPatchSizeTest has been previously run. 

for iProj = 1:numel(projList)
    load([projList{iProj} filesep 'ANALYSIS' filesep 'movieData.mat']);
    movieData= MD;
    threshDir = [movieData.outputDirectory_ filesep ....
        'neurite_body_masks' filesep 'LocalThreshPatchSizeTest'];
    files = searchFiles('.png',[],threshDir,0,'all',1);
    
    % cat
    imgs = cellfun(@(x) imread(x),files,'uniformoutput',0);
    imgLarge = horzcat(imgs{:});
    fsFigure(.75);
    imshow(imgLarge);
    
    h = msgbox('Examine Local Threshold');
    uiwait(h)
    
    %
    prompt = {'Enter Local Threshold Patch Size (In Pixels)'};
    def = {'40'};
    
    patchSize = inputdlg(prompt,'Local Threshold Input',1,def);
    save([threshDir filesep 'manualPatchSizeSelect.mat'],'patchSize');
    close all
    
end
end

