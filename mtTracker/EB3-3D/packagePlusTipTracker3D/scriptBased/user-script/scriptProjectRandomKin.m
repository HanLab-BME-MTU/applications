%%
MLPath='/home2/proudot/Danuser_lab/externBetzig/analysis/proudot/anaProject/sphericalProjection/';
movieListFileNames={'prometaphase/analysis/ML1min.mat','prometaphase/analysis/ML4min.mat'};
%% movieListFileNames={'prometaphase/analysis/fourProm.mat'};
% Build the array of MovieList (automatic)
aMovieListArray=cellfun(@(x) MovieList.loadMatFile([MLPath filesep x]), movieListFileNames,'unif',0);
aMovieListArray=[aMovieListArray{:}];
%%

for mlIdx=1:length(aMovieListArray)
    ML=aMovieListArray(mlIdx);
    bundledTracks=cell(1,length(ML.movieDataFile_));
    bundledTracksRandom=cell(1,length(ML.movieDataFile_));
    [handles,~,hf]=setupFigure(length(ML.movieDataFile_),1,length(ML.movieDataFile_));
    set(hf,'name',movieListFileNames{mlIdx},'numbertitle','off')

    %% loop over each different cell in each condition
    for mIdx=1:length(ML.movieDataFile_)
        MD=MovieData.loadMatFile(ML.movieDataFile_{mIdx});
        %bundleStatisticsWholeSpindleTestingInlier(MD);
        outputDirBundle=[MD.outputDirectory_ filesep 'Kin' filesep 'bundles' ];

        name='inliers';
        tmp=load([outputDirBundle filesep name filesep  'kin-MT-bundle.mat']);
        kinTracks=tmp.kinTracks;
        bundledTracks{mIdx}=kinTracks;
        
        name='inliersRandom';
        tmp=load([outputDirBundle filesep name filesep  'kin-MT-bundle.mat']);
        kinTracksRandom=tmp.kinTracks;
        bundledTracksRandom{mIdx}=kinTracksRandom;
    
        
        %% Number of fiber at each given frame. 
        kinTracksCell={bundledTracks{mIdx},bundledTracksRandom{mIdx}};
        fiberCount=zeros(length(kinTracksCell),kinTracksCell{1}.numTimePoints);
        kinetochoreCount=zeros(length(kinTracksCell),kinTracksCell{1}.numTimePoints);
        
        for tIdx=1:length(kinTracksCell)
            kinTracks=kinTracksCell{tIdx};
            for kIdx=1:length(kinTracks)
                fiberCount(tIdx,kinTracks(kIdx).f)=fiberCount(tIdx,kinTracks(kIdx).f)+double(length(kinTracks(kIdx).fiber));
                kinetochoreCount(tIdx,kinTracks(kIdx).f)=kinetochoreCount(tIdx,kinTracks(kIdx).f)+1;
            end
        end
        plot(handles(mIdx),linspace(0,kinTracksCell{1}.numTimePoints*MD.timeInterval_,kinTracksCell{1}.numTimePoints), fiberCount./kinetochoreCount);
        xlabel(handles(mIdx),'Frame count');
        ylabel(handles(mIdx),'avg MT/bundle');
        ylim(handles(mIdx),[10 35]);
        legend(handles(mIdx),{'Kin','Random'});
    end
end


