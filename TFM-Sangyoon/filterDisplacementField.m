function [displField] = filterDisplacementField( orgDisplacementField, mask)
%cropMovieTFMPackage crops movie images and displacement field for finer
%force reconstruction.
% input:    pathForTheMovieDataFile:    path to the movieData file
% output:   cropped displacement stored in pathForTheMovieDataFile/correctedDisplacementField
%           original displacement stored in pathForTheMovieDataFile/correctedDisplacementField
%           cropped images stored at a folder where their orinials were
%           original channels are moved to one folder inside its own folder
% Sangyoon Han April 2013
    % Get whole frame number
    nFrames = length(orgDisplacementField);
    displField(nFrames)=struct('pos','','vec','');
    for k=1:nFrames
        % filtering out the flow vectors outside the mask
        idx = true(1,numel(orgDisplacementField(k).pos(:,1)))'; %index for filtering

        outsideIdx = arrayfun(@(i,j) maskVectors(i,j,mask),orgDisplacementField(k).pos(:,1),orgDisplacementField(k).pos(:,2));
        idx = idx & outsideIdx;
        displField(k).pos = orgDisplacementField(k).pos(idx,:);
        displField(k).vec = orgDisplacementField(k).vec(idx,:);
    end
    
end