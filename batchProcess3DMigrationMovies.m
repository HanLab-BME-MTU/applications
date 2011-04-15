function [movieArray,errMess] = batchProcess3DMigrationMovies(movieArray,varargin)
%BATCHPROCESS3DMIGRATIONMOVIES runs all analysis steps on every input 3D migration movie
% 
%   errMess = batchProcess3DMigrationMovies(movieArray)
%   errMess = batchProcess3DMigrationMovies(movieArray,'OptionName1',optionValue1,...)
% 
% This function runs process3DMigrationMovie on every input 3d movie, in
% parallel. Each movie's results for each step are saved to its analysis
% directory. Movies which cause errors will be skipped, and the error
% messages returned. The same parameters will be used for every movie, and
% each movie must have the same number and order of channels.
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
%       ('NumParallel' -> positive integer scalar) Specifies the number of
%       movies to run simultaneously. If the job is run locally, the max is
%       8. If run on orchestra, it's ... a lot???
%       Default is 8.
%
%       *NOTE:* Additionally, any options supported by
%       process3DMigrationMovie can be input, and will be passed to this
%       function when it is called.
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
p.addParamValue('NumParallel',8,@(x)(numel(x) == 1));

p.parse(movieArray,varargin{:});
p = p.Results;


%% -----------------------------Init --------------------------------- %%


%Run sanity check on all movie datas, and relocate if necessary.




nMovies = numel(movieArray);

%Check if we have LSF available, and try to set it up:
try
    disp('Attempting to configure LSF....')
    lsf = findResource('scheduler','type','lsf');
    lsf.ClusterMatlabRoot = '/opt/matlab';
    lsf.DataLocation = '/groups/lccb-nih';
    job = createParallelJob(lsf);
    
    
    
    doLocal = false;
catch em
    %Otherwise, just run it locally.
    disp(['Running processing locally - LSF scheduler could not be configured  : ' em.message]);
    doLocal = true;
    if p.NumParallel > 8
        disp('Maximum of 8 workers available locally - setting NumParallel to 8!')
    end
end

if doLocal
    disp('Setting up local parallel job...')
    %Set up the local workers
    poolSize = matlabpool('size');
    if poolSize ~= p.NumParallel    
        if poolSize > 0
            matlabpool('close');
        end
        matlabpool('open',p.NumParallel);
    end
else    
    
    
    
    
    
    
end
   


