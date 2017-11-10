function [ out ] = saveStats( ML )
%saveStats saves statitics on lamins to disk

ML = lamins.util.toMovieList(ML);
out.ML = ML;
out.area = lamins.functions.getFaceStatistics(ML,'Area',1000);
out.eccentricity = lamins.functions.getFaceStatistics(ML,'Eccentricity',1);
out.perimeter = lamins.functions.getFaceStatistics(ML,'Perimeter',sqrt(1000));
out.edgeLength = lamins.analysis.getEdgeStatistics(ML,'Area',sqrt(1000));
out.edgesPerVertex = lamins.analysis.getEdgeStatistics(ML,'EdgesPerVertex',1);
out.edgesPerFace = lamins.functions.getFaceStatistics(ML,'EdgesPerFace',1);
out.nucleusCircularity = lamins.analysis.getCircularity(ML);
save([ML.movieListPath_ filesep 'stats_20150622.mat'],'-struct','out');


end

