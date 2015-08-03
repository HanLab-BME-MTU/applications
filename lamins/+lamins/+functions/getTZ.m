function [ tzs ] = getTZ( ML, idx, filename )
% Obtains the tz variable from the saved output of an analyzeLamins script

if(nargin < 2)
    idx = 1;
end
if(nargin < 3)
    filename = 'skeletons_2015_06_10.mat';
end

if(isempty(ML.movies_))
    ML.sanityCheck();
end

directories = cellfun(@(MD) MD.outputDirectory_,ML.movies_,'Unif',false)';
data = cellfun(@(D) load([D filesep filename],'tz'),directories,'Unif',false);
tzs = cellfun(@(data) data.tz(idx),data,'ErrorHandler',@(s,in) Inf);

end

