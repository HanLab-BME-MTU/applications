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
    disp('Start filetering field with roiMask...'); tic
    nFrames = length(orgDisplacementField);
    displField(nFrames)=struct('pos','','vec','');
    outsideIdx = arrayfun(@(y) maskVectors(y.pos(:,1),y.pos(:,2),mask),orgDisplacementField,'UniformOutput',false);
    for k=1:nFrames
        displField(k).pos = orgDisplacementField(k).pos(outsideIdx{k},:);
        displField(k).vec = orgDisplacementField(k).vec(outsideIdx{k},:);
    end
    toc
    disp('Done.')
end
% old code - very inefficient... commented out on 5/9/16
%     for k=1:nFrames
%         % filtering out the flow vectors outside the mask
%         idx = true(1,numel(orgDisplacementField(k).pos(:,1)))'; %index for filtering
% 
% %         outsideIdx = arrayfun(@(i,j) maskVectors(i,j,mask),orgDisplacementField(k).pos(:,1),orgDisplacementField(k).pos(:,2));
%         outsideIdx = arrayfun(@(y) maskVectors(y(:,1),y(:,2),mask),orgDisplacementField(k).pos);
%         idx = idx & outsideIdx;
%         displField(k).pos = orgDisplacementField(k).pos(idx,:);
%         displField(k).vec = orgDisplacementField(k).vec(idx,:);
%     end
