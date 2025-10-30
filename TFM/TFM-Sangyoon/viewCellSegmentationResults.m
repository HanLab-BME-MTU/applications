% viewAdhesionSegmentationImagesFromAdhAnalBatchResults go through each
% movie and open the segmented images

%% load ML
[pathAnalysisAll, MLNames, groupNames, usedSelectedFoldersMat,...
    specificName,~,MLdirect]=chooseSelectedFolders;

%% Load movieLists for each condition
numConditions = numel(pathAnalysisAll);
for k=1:numConditions
    MLAll(k) = MovieList.load([pathAnalysisAll{k} filesep MLNames{k}]);
end
%% Each ML
for ii=1:numConditions
    N(ii) = numel(MLAll(ii).movies_);
end
jformat = ['%.' '3' 'd'];
k=0;
%% Go through each movie
k=k+1;
iML=0;
curN = 0;
% nP =0;
for p=1:k
%     nP=nP+1;
    if p>(curN+N(iML+1))
        iML=iML+1;
        curN = curN+N(iML);
    end
end
iMD = k-curN;

curML = MLAll(iML+1);
curMD = curML.movies_{iMD};

maskProc = curMD.getProcess(curMD.getProcessIndex('MaskRefinementProcess'));
iChan  = find(maskProc.checkChannelOutput);
curImg = curMD.channels_(iChan).loadImage(1);
figure, imshow(curImg,[]), hold on
mask = maskProc.loadChannelOutput(iChan,1);
[B,~,nBD]  = bwboundaries(mask,'noholes');
for kk=1:nBD
    boundary = B{kk};
    plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 2) % cell boundary
end


disp('Current movie is: ')
disp(curMD.outputDirectory_)
disp('Movie name is: ')
[~,curFName]=fileparts(curMD.outputDirectory_);
disp(curFName)

