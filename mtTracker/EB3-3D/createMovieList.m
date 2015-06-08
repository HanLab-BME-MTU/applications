function movieListArray=createMovieList(dataRootPath,varargin)
% Create a MovieList for each set of movies associated to a or condition. 
% Contract: 
%  - Element in <condNames> must correpond to folders name in which cell movies are gathered
%  - Each movies must be organized following:
%      / <dataRootPath> / <condNames{i}> / <cellFolder> / <channelFolder> / time-point-0001.tif

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataRootPath', @ischar);
ip.addParamValue('timeInterval',1, @isnumeric);
ip.addParamValue('lateralPixelSize',200, @isnumeric);
ip.addParamValue('axialPixelSize',1, @isnumeric);
ip.addParamValue('condNames',{}, @iscell);
ip.addParamValue('resPath',dataRootPath, @ischar);
ip.addParamValue('pixelSize',[1 1], @ischar);
ip.parse(dataRootPath,varargin{:});

resPath=ip.Results.resPath;
condNames=ip.Results.condNames;
movieListArray=MovieList.empty(0,length(condNames));
if(isempty(condNames))
    condFolder=dir([dataRootPath]);
    condFolder = condFolder(arrayfun(@(x) (x.isdir)&&(x.name(1)~= '.'), condFolder) );
    condNames={condFolder.name};
end

for i=1:length(condNames)
    mkdir([resPath filesep condNames{i}]);
    cellFolder=dir([dataRootPath filesep condNames{i} filesep]);
    cellFolder = cellFolder(arrayfun(@(x) (x.name(1)~='.')&(x.isdir), cellFolder));
    
    ML=cell(1,length(cellFolder));
    for k=1:length(cellFolder)
        if (cellFolder(k).isdir)
            channelFolder=dir([dataRootPath filesep condNames{i} filesep cellFolder(k).name]);
            channelFolder = channelFolder(arrayfun(@(x) x.name(1)~= '.', channelFolder) );
            channelCell=Channel.empty(0,length(channelFolder));
            for j=1:length(channelFolder)
                try 
                    channelCell(j)=Channel([dataRootPath filesep condNames{i} filesep cellFolder(k).name filesep channelFolder(j).name]);
                    channelCell(j).getImageFileNames;
                catch 
                    channelCell(j)=Channel();
                end 
            end
            channelCell=channelCell(arrayfun(@(x) ~isempty(x.channelPath_),channelCell));
            outputDir=[resPath filesep condNames{i} filesep cellFolder(k).name filesep 'analysis'];
            mkdir(outputDir);
             try
                disp([dataRootPath filesep cellFolder(k).name]);
                MD=MovieData(channelCell,outputDir, ...
                    'pixelSize_',ip.Results.lateralPixelSize,'pixelSizeZ_',ip.Results.axialPixelSize, ...
                    'timeInterval_',ip.Results.timeInterval, ...
                    'movieDataPath_' , outputDir, ...
                    'movieDataFileName_', 'movieData.mat');
                MD.sanityCheck();
                ML{k}=MD;
             catch
                 error(['Movie ' cellFolder(k).name ' failed to build']);
             end
        end
    end
    try
    MList=MovieList([ML{:}],[resPath filesep condNames{i}],'movieListPath_',[resPath filesep condNames{i}],'movieListFileName_','movieList.mat');
    MList.sanityCheck();
    movieListArray(i)=MList;
    catch
        error(['Condition ' condNames{i} ' failed to build']);
    end 
end




