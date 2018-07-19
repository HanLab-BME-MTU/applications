function T=collectDetectionData(movieList)
% Summarize detection information in a single CSV file
% P. Roudot 2016
ip=inputParser();
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('movieList', @(x) isa(x,'MovieList'));
ip.parse(movieList);
p=ip.Results;

movieList.sanityCheck();

T=table();
for mIdx=1:length(movieList.movieDataFile_)
    %%
    MD=movieList.getMovie(mIdx);
    %identify the last detection process
    detectionProcesses=cellfun(@(p) isa(p,'DetectionProcess'),MD.processes_);
    pIdx=find(detectionProcesses,1,'last');
    
    if(~isempty(pIdx))
        pr=MD.getProcess(pIdx);
        movieTable=table();
        for oIdx=1:length(pr.outFilePaths_)
            timeInter=1;
            if(~isempty(MD.timeInterval_)) timeInter=MD.timeInterval; end
            movieTable=[movieTable;detection2table(pr.loadChannelOutput(oIdx),timeInter,oIdx)];
        end
        [~,f,~]=fileparts(MD.movieDataFileName_);
        name=cell(height(movieTable),1);
        [name{:}]=deal(f);
        movieTable.movieId=mIdx*ones(height(movieTable),1);
        movieTable.movieName=name;
        T=[T;movieTable];
    else
        warning(['No detection process found in Movie Data ' MD.movieDataFileName_ ]);
    end
    
    
end
%reorder table
T=[T(:,end-1:end) T(:,1:end-2)];
[~,f,~]=fileparts(movieList.movieListFileName_);
writetable(T,[movieList.outputDirectory_ filesep f '_detections.csv']);

