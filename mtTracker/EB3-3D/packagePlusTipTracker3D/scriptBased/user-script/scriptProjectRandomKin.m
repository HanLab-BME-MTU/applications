%%
MLPath='/home2/proudot/Danuser_lab/externBetzig/analysis/proudot/anaProject/sphericalProjection/';
movieListFileNames={'prometaphase/analysis/fourProm.mat'};
%% movieListFileNames={'prometaphase/analysis/fourProm.mat'};
% Build the array of MovieList (automatic)
aMovieListArray=cellfun(@(x) MovieList.loadMatFile([MLPath filesep x]), movieListFileNames,'unif',0);
aMovieListArray=aMovieListArray{:};
%%

for k=1:length(aMovieListArray)
    ML=aMovieListArray(k);
    bundledTracks=cell(1,length(ML.movieDataFile_));
    bundledTracksRandom=cell(1,length(ML.movieDataFile_));

    %% loop over each different cell in each condition
    for i=1:length(ML.movieDataFile_)
        MD=MovieData.loadMatFile(ML.movieDataFile_{i});
        %bundleStatisticsWholeSpindleTestingInlier(MD);
        outputDirBundle=[MD.outputDirectory_ filesep 'Kin' filesep 'bundles' ];

        name='inliers';
        tmp=load([outputDirBundle filesep name filesep  'kin-MT-bundle.mat']);
        kinTracks=tmp.kinTracks;
        bundledTracks{i}=kinTracks;
        name='inliersRandom';
        tmp=load([outputDirBundle filesep name filesep  'kin-MT-bundle.mat');
        kinTracks=tmp.kinTracks;
        bundledTracksRandom{i}=kinTracks;
        [handles,hFig]=bundleStatistics(MD,'kinBundle',{bundledTracks{i},bundledTracksRandom{i}},'kinBundleName',{'Kinetochore','Randomized Kinetochore'});
%         title(hFig,
    end
end


