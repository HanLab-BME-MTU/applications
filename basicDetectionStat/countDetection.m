function countDetection(movieListCell)
ip=inputParser();
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('movieListCell', @iscell);
ip.parse(movieListCell);
p=ip.Results;

maxChannel=5;
cellcount=cell(maxChannel,length(movieListCell));
for cIdx=1:length(movieListCell)
    movieListCell{cIdx}.sanityCheck();
end
for cIdx=1:length(movieListCell)
    ML=movieListCell{cIdx};
   %     movieCount=nan(1,length(ML.movieDataFile_));
    disp(['Movie List: ' ML.getFilename]);
    detectedMovieCount=1;
    for mIdx=1:length(ML.movieDataFile_)
        %%
            MD=ML.getMovie(mIdx);
            %identify the last detection process
            detectionProcesses=cellfun(@(p) isa(p,'DetectionProcess'),MD.processes_);
            pIdx=find(detectionProcesses,1,'last');
            
            if(~isempty(pIdx))
                fprintf([ MD.movieDataFileName_]);
                pr=MD.getProcess(pIdx);
                for oIdx=1:length(pr.outFilePaths_)
                    cellcount{oIdx,cIdx}(detectedMovieCount)=length(pr.loadChannelOutput(oIdx).xCoord);                  
                    fprintf([' channel ' num2str(oIdx) ': ' num2str(cellcount{oIdx,cIdx}(detectedMovieCount))]);
                end
                detectedMovieCount=detectedMovieCount+1;
                fprintf('\n');
            else
                warning(['No detection process found in Movie Data ' MD.movieDataFileName_ ]);
            end
                
            
    end
    fprintf('\n');
%     cellcount{cIdx}=movieCount;
end

for chIdx=1:maxChannel
    if(~isempty(cell2mat(cellcount(chIdx,:))))
        figure();
        grouping=cellfun(@(c,i) i*ones(size(c)),cellcount(chIdx,:),num2cell(1:length(cellcount(chIdx,:))),'unif',0);
        boxplot(cell2mat(cellcount(chIdx,:)),cell2mat(grouping));
        ylim([0,1.1*max(cell2mat(cellcount(chIdx,:)))]);
        ylabel('Detection Count')
        title(['Channel ' num2str(chIdx)]);
    end
end
