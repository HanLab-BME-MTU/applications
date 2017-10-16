function MD=createMD(channelPaths,pixelSize,pixelSizeZ,timeInterval)
chs=cellfun(@(c) Channel(c),channelPaths,'unif',0);
outputFolder=fileparts(channelPaths{1});
mkdirRobust([outputFolder filesep 'analysis']);
MD=MovieData([chs{:}],[outputFolder filesep 'analysis'],'movieDataFileName_','movieData.mat','movieDataPath_',[outputFolder filesep 'analysis'], ...
	    'pixelSize_',pixelSize,'pixelSizeZ_',pixelSizeZ,'timeInterval_',timeInterval);
MD.save();