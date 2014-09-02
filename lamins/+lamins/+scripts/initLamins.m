import lamins.classes.*;
import lamins.functions.*;

currentDir = pwd;
cd('~/matlab/lamins');
ML = MovieList.load('project/140515/original reconstructed images/movieList.mat');
Lamins(1) = LaminsData(ML.getMovie(1));
Lamins(2) = LaminsData(ML.getMovie(2));
cd(currentDir);

%% Generate image classes (quick)
I{1} = squeeze(Lamins(1).getImages);
I{2} = squeeze(Lamins(2).getImages);
