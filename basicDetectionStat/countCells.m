function T=countCells(varargin)
% Count the total number of detection in the movie per-channel.
% P. Roudot 2016
ip=inputParser();
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addOptional('movieListDetectionFile',[],@ischar);
ip.parse(varargin{:});
p=ip.Results;

if(isempty(p.movieListDetectionFile))
    [basename,folder] = uigetfile({'*.csv;*.mat'});
    p.movieListDetectionFile=fullfile(folder,basename);
end
[filePath,filename,ext]=fileparts(p.movieListDetectionFile);
if(strcmp(ext,'.csv'))
    T=readtable(p.movieListDetectionFile);
elseif(strcmp(ext,'.mat'))
    disp('Collecting movie detection.')
    movieList=MovieList.load(p.movieListDetectionFile);
    T=collectDetectionData(movieList);
    filePath=[movieList.outputDirectory_];
    filename=[filename '_detections'];
end
    
disp('Processing detection.')
%%
T.count=T.C*0;
% For each movie ID in the table
for movieId=unique(T.movieId)'  
    movieRow=(T.movieId==movieId);
    % For each channel in the table
    for cIdx=unique(T.C)'
       % identify the detections belonging to those movies and channel 
        selectedRow=movieRow & (T.C==cIdx);
       % compute their sum and selectively add them to the count row
       T.count=T.count + selectedRow*sum(selectedRow);
    end
end
% Keep only the value we are interested in
T=T(:,{'movieId','movieName','C','count'});
% Remove redundancy in the table
T=unique(T);
writetable(T,[filePath filesep filename '_counts.csv']);



