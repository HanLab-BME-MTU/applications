tmp=load('/work/gdanuser/proudot/project/EB3-3D-track/packaging/alpha/plusTipTracker3D-alpha-1/data/A1_HeLa_Cells_EB1_and_H2B/xpMovieListCollection.mat');
ML=tmp.movieListCell{1}
crop3D(ML.getMovies{1},[129,163,19,367,442,140]);