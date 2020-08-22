% createSingleFrameMDsFromMultiFrameMDs.m
% This script sample only the first frame from the movies and a new
% movieData with the same MD parameters.
% Sangyoon Han Aug 2020
%% Read selectedFolders.mat
[pathAnalysisAll, MLNames, groupNames,usedSelectedFoldersMat] = chooseSelectedFolders;
%% Load movieLists for each condition
numConditions = numel(pathAnalysisAll);
for k=1:numConditions
    MLAll(k) = MovieList.load([pathAnalysisAll{k} filesep MLNames{k}]);
end
%% Making single frame MD
for ii=1:numConditions
    curML=MLAll(ii);
    curMovies = curML.movies_;
    N = numel(curMovies);
    
    for k=1:N % per movie
        % Load the movie
        cellMD = curMovies{k};
        % Setup storage location
        pathOneFrame = [cellMD.outputDirectory_ 'OneFrame'];
        mkdir(pathOneFrame)
        
        % Sample the first image per channels
        for iiChan=1:numel(cellMD.channels_)
            curChan=cellMD.channels_(iiChan);
            curFirstImg = curChan.loadImage(1); 
            %Create a new folder with one frame
            path3 = [pathOneFrame filesep num2str(curChan.excitationWavelength_)];
            mkdir(path3)
            firstImgPath = [path3 filesep num2str(curChan.excitationWavelength_) '.tif'];
            imwrite(uint16(curFirstImg),firstImgPath,'Compression','none')
            % Register as channels
            curChanNew(iiChan)=Channel(path3);
            curChanNew(iiChan).excitationWavelength_ = curChan.excitationWavelength_;
            curChanNew(iiChan).emissionWavelength_ = curChan.emissionWavelength_;
        end
        % Make into an MD
        cellMD2(k) = MovieData(curChanNew,pathOneFrame);
        % Copy MD metadata
        cellMD2(k).setPath(pathOneFrame);
        cellMD2(k).setFilename('movieData.mat');
        cellMD2(k).numAperture_=cellMD.numAperture_;
        cellMD2(k).camBitdepth_=cellMD.camBitdepth_;
        cellMD2(k).timeInterval_ = cellMD.timeInterval_;
        cellMD2(k).pixelSize_= cellMD.pixelSize_; % 60x x 1.8x (new objective config.)
        cellMD2(k).sanityCheck;
        cellMD2(k).save        
    end
    ML = MovieList(cellMD2,pathAnalysisAll{ii});
    ML.setPath(pathAnalysisAll{ii});
    MLNames{ii} = 'movieListOneFrame.mat';
    ML.setFilename(MLNames{ii});
    ML.sanityCheck;
    ML.save
    
    [~,groupNames{ii}] = fileparts(pathAnalysisAll{ii});
    groupNames{ii}=groupNames{ii}(end-10:end);
end
rootAnalysis=pathAnalysisAll{1};
specificName = strjoin(groupNames);
save([rootAnalysis filesep 'selectedFoldersOneFrame' specificName '.mat'], 'rootAnalysis','pathAnalysisAll','MLNames','groupNames')
