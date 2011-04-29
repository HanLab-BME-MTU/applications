function runInfo = batchProcess3DMigrationMovies(movieArray,varargin)
%BATCHPROCESS3DMIGRATIONMOVIES runs all analysis steps on every input 3D migration movie
% 
%   runInfo = batchProcess3DMigrationMovies(movieArray)
%   runInfo = batchProcess3DMigrationMovies(movieArray,'OptionName1',optionValue1,...)
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
%   runInfo - An array the same size as movieArray containing information
%   on the processing of each movie, including any error messages, the time
%   required to complete, parameters used etc.
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
ip = inputParser;
ip.FunctionName = mfilename;
ip.KeepUnmatched = true; %Keep extra parameters for passing to processing function

ip.addRequired('movieArray',@(x)(isa(x,'MovieData3D')));
ip.addParamValue('NumParallel',8,@(x)(numel(x) == 1));

ip.parse(movieArray,varargin{:});
p = ip.Results;%Parameters for this function
fp = ip.Unmatched;%Parameters to pass to processing function


%% -----------------------------Init --------------------------------- %%

nMovies = numel(movieArray);

%Set up general default parameters used on all movies.
if isfield(fp,'BatchMode') && fp.BatchMode == false
    %if the user tried to turn it off, warn them that this isn't supported
    warning('Migration3D:batchProcess:batchMode','Batch mode must be enabled to run processing in parallel!')
end
fp.BatchMode = true;
%replicate parameters for passing to the processing function for each movie
fp = repmat(fp,nMovies,1);


%Run sanity check on all movie datas, and relocate if necessary.
%NOT DOING YET. NEED TO INPUT MOVIELIST instead AND THEN LOAD MOVIE array - TEMP


%Check if we have LSF available, and try to set it up:
try
    disp('Attempting to configure LSF....')
    lsf = findResource('scheduler','type','lsf');
    lsf.ClusterMatlabRoot = '/opt/matlab';
    lsf.DataLocation = '/groups/lccb-nih';
    job = createJob(lsf);        
    
    %TEMP, under construction!
    
    doLocal = false;
catch em
    %Otherwise, just run it locally.
    disp(['Running processing locally - LSF scheduler could not be configured  : ' em.message]);
    doLocal = true;
    if p.NumParallel > 8
        disp('Maximum of 8 workers available locally - setting NumParallel to 8!')
        p.NumParallel = 8;
    end  
end

%% -------------------------- Processing --------------------------- %%


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
    
    runInfo = cell(nMovies,1);
    
    parfor iMov = 1:nMovies
                
        try            
            runInfo{iMov} = process3DMigrationMovie(movieArray(iMov),fp(iMov));                         
        catch em
            disp(['Error processing movie ' num2str(iMov) ' : ' em.message])
            runInfo{iMov} = em;
        end
            
    end
    
else    
    
    
    disp('not set up yet!')
    
    
    
    
end
   


