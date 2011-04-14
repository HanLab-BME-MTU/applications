function [movieArray,errMess] = batchProcess3DMigrationMovies(movieArray,varargin)
%BATCHPROCESS3DMIGRATIONMOVIES runs all analysis steps on every input 3D migration movie
% 
%   errMess = batchProcess3DMigrationMovies(movieArray)
%   errMess = batchProcess3DMigrationMovies(movieArray,'OptionName1',optionValue1,...)
% 
% This function runs all required processing/analysis steps on every input
% 3D migration movie: segmentation, skeletonization, pruning, etc. etc. It
% will analyze multiple movies in parallel. Each movie's results for each
% step are saved to its analysis directory. Movies which cause errors will
% be skipped, and the error messages returned.
% 
% 
% Input:
% 
%   movieArray - an array of MovieData3D objects. Every input movie must
%   have the same number and order of channels. Can be created with
%   setupMovieArray.m
%
%   'OptionName',optionValue - A string with an option name followed by the
%   value for that option.
% 
%   Possible Option Names:
%
%       ('OptionName' -> possible values)
%
%       ('SegChannelIndex -> positive integer scalar) The index of the
%       channel to use for segmentation.
%       Default is 1.
%       
%       ('RunOrchestra' -> true/false) If true, the job will be submitted
%       to Orchestra for processing there. If false, it will be run
%       locally.
%       Default is false.
% 
% Output:
%
%   errMess - An array the same size as movieArray containing any error
%   messages which were generated during the processing of the movies.
%
%   Additionally, the results of each processing/analysis step will be
%   saved to each movie's analysis directory, and logged in the processes
%   array of the MovieData3D object.
%
% Hunter Elliott
% 4/14/2011
%

%% --------------------------- Input --------------------------------- %%

%Parse all the inputs
p = inputParser;
p.FunctionName = mfilename;

p.addRequired('movieArray',@(x)(isa(x,'MovieData3D')));
p.addParamValue('SegChannelIndex',1,@(x)(numel(x) == 1 && isposint(x)));
p.addParamValue('RunOrchestra',false,@(x)(numel(x) == 1));

p.parse(movieArray,varargin{:});


%% -----------------------------Init --------------------------------- %%


poolSize = matlabpool('size');
if poolSize ~= p.NumParallel    
    if poolSize > 0
        matlabpool('close');
    end
    matlabpool('open',p.NumParallel);
end



