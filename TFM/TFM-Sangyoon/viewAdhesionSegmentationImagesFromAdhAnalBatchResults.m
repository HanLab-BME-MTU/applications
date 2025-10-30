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
for p=1:k
    if p>(curN+N(iML+1))
        iML=iML+1;
        curN = curN+N(iML);
    end
end
iMD = k-curN;

curML = MLAll(iML+1);
curMD = curML.movies_{iMD};

outputFilePath = [curMD.outputDirectory_ filesep 'Adhesion Quantification'];
imgPath = [outputFilePath filesep 'imgs'];
figPath = [imgPath filesep 'figs'];
%will open only first one
fig=openfig(strcat(figPath,'/imgNAFA',num2str(1,jformat)));

disp('Current movie is: ')
disp(curMD.outputDirectory_)
disp('Movie name is: ')
[~,curFName]=fileparts(curMD.outputDirectory_);
disp(curFName)
%% save to designated place
curAx = fig.Children;
%Add a scale bar
line([curAx.XLim(1)+5 curAx.XLim(1)+5+5000/curMD.pixelSize_],...
    [curAx.YLim(1)+10 curAx.YLim(1)+10],'linewidth',3,'color','w')
saveDir = '/storage/network/Collaboration/AnneCress_UArizona/MultiImageAnalysis/repImgs';
if ~exist(saveDir,'dir')
    mkdir(saveDir)
end
savefig(fig,[saveDir filesep curFName 'seg.fig'])
