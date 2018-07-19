function [meanStdFeatNum,idx2rem,p]=plusTipCheckFeatureDensity(projList,prctl2rem)
% checkFeatureDensity checks whether feature number changes significantly over the course of a movie
%
% SYNOPSIS: [meanStdFeatNum,idx2rem,p]=plusTipCheckFeatureDensity(projList,prctl2rem)
%
% INPUT : projList  : list of projects, output of getProj
%         prctl2rem : percentile of std above which movies should be
%                     removed
%         p         : value corresponding to the percentile
% OUTPUT: meanStdFeatNum: nProj x 2 matrix containing the mean and standard
%                         deviation for the number of detected features
%                         over all frames of the movie.  A high standard
%                         deviation may indicate significant photobleaching
%                         or focus drift.
%         idx2rem       : indices of movies where the std is in the top
%                         prctl2rem-percentile
%         histogram - of std of normalized feature number across frames,
%         with the red line representing the percentile given



nMovies=size(projList,1);
meanStdFeatNum=zeros(nMovies,2);

for i=1:nMovies
    % load detection file movieInfo.mat and convert to cell array
    temp=load([formatPath(projList(i,1).anDir) filesep 'feat' filesep 'movieInfo.mat']);
    movieInfo=temp.movieInfo;
    movInfoCell=struct2cell(movieInfo);
    nFrames=length(movieInfo);
    featNumPerFrame=nan(nFrames,1);
    idx=find(~cellfun(@isempty, movInfoCell(1,:)));
    % get number of features in each frame
    featNumPerFrame(idx)=cell2mat(cellfun(@(x) length(x(:,1)),movInfoCell(1,idx),'UniformOutput',0))';
    % normalize by the mean feature number for this movie
    %featNumPerFrame=featNumPerFrame./nanmean(featNumPerFrame);
    % store data
    meanStdFeatNum(i,:)=[nanmean(featNumPerFrame) nanstd(featNumPerFrame)];
end

p=prctile(meanStdFeatNum(:,2),prctl2rem);
idx2rem=find(meanStdFeatNum(:,2)>p);

figure
[n,xout] = hist(meanStdFeatNum(:,2),20);
hist(meanStdFeatNum(:,2),20)
hold on
plot([p;p],[0;max(n)],'r')
xlabel('std of normalized feature number across frames')
ylabel('number of movies')