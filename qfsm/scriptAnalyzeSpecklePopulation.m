% Load MD
laptopPath = '';
serverPath='/storage/network/ActinSpeckleMicroscopy/CodeDevelopment/Speckle11_TalinDCTS_Cell_02/movieData.mat';
MD = MovieData.load(serverPath);

analyzeSpecklePopulation(MD)