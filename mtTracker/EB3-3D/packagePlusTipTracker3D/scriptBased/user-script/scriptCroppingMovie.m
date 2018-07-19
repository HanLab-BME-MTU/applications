root='/project/cellbiology/gdanuser/december/philippe/externBetzig/analysis/proudot/anaProject/nodeEstimationBug/';

MD=MovieData.load([root filesep 'cell1_6/deskew/analysis/movieData.mat']);
crop3D(MD,[50,1,1,MD.imSize_(1),MD.imSize_(2),MD.zSize_]);

MD=MovieData.load([root filesep 'cell1_9/deskew/analysis/movieData.mat']);
crop3D(MD,[1,1,1,MD.imSize_(1)-50,MD.imSize_(2),MD.zSize_]);

croppedMovieList=indexLatticeData([root filesep 'cell1_*/deskew/analysis/cropped/ch{ch}/*tif'],root);