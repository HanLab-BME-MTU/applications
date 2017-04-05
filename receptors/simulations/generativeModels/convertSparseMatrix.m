for i=1:length(compTracks)
    compTracks(i).tracksFeatIndxCG = full(compTracks(i).tracksFeatIndxCG);
    compTracks(i).tracksCoordAmpCG = full(compTracks(i).tracksCoordAmpCG); 
    compTracks(i).tracksCoordAmpCG(compTracks(i).tracksCoordAmpCG==0) = NaN;
end