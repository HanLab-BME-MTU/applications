function [meanStdFeatNum]=plusTipCheckFeatureDensity(projList)
% checkFeatureDensity checks whether feature number changes significantly over the course of a movie

% INPUT : projList: list of projects, output of getProj
% OUTPUT: meanStdFeatNum: nProj x 2 matrix containing the mean and standard
%                         deviation for the number of detected features
%                         over all frames of the movie.  A high standard
%                         deviation may indicate significant photbleaching
%                         or focus drift.



nMovies=size(projList,1);
meanStdFeatNum=zeros(nMovies,2);

for i=1:nMovies
    % load detection file movieInfo.mat and convert to cell array
    temp=load([formatPath(projList(i,1).anDir) filesep 'feat' filesep 'movieInfo.mat']);
    movieInfo=temp.movieInfo;
    movInfoCell=struct2cell(movieInfo);
    % get number of features in each frame
    featNumPerFrame=cell2mat(cellfun(@(x) length(x(:,1)),movInfoCell(1,:),'UniformOutput',0))';
    % normalize by the mean feature movie for this frame
    featNumPerFrame=featNumPerFrame./mean(featNumPerFrame);
    % store data
    meanStdFeatNum(i,:)=[mean(featNumPerFrame) std(featNumPerFrame)];
end

figure
hist(meanStdFeatNum(:,2),20)
xlabel('std of normalized feature number across frames')
ylabel('number of movies')
