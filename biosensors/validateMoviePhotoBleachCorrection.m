function movieData = validateMoviePhotoBleachCorrection(movieData,varargin)

% movieData = validateMoviePhotoBleachCorrection(movieData)
% 
% movieData = validateMoviePhotoBleachCorrection(movieData,'OptionName',optionValue,...)
% 
% This function checks the quality of the fit used in the photobleach
% correction for the input movie. The results of this validation are stored
% in the movieData.
% 
% Input: 
% 
%   movieData - The structure describing the movie, as created with
%   setupMovieData.m
% 
%   'OptionName',optionValue - A string with an option name followed by the
%   value for that option.
% 
%   Possible Option Names:
%
%       ('OptionName' -> possible values)
%
%       ('ChannelIndex -> Positive integer scalar or vector) The integer
%       index of the channel to validate the photobleach
%       correction on. This index corresponds to the channel directory
%       location in the cell array movieData.channelDirectory.
%       If not input, the user will be asked to select a channel.
%
%       ('PortTest' -> True/False) If true, the residuals of the photobleach
%       correction fit will be tested for temporal correlation using the
%       portmanteau test. 
%       Default is true.
%
%       ('NormTest' -> True/False) If true, the residuals of the
%       photobleach correction fit will be tested for normality using the
%       Kolmogorov-Smirnov test.
%       Default is true.
%
%       ('Alpha' -> Scalar between 0 and 1) The significance level to
%       perform all tests at.
%       Default is .05
%
%
%  NOTE: If multiple tests are specified, the fit must pass ALL test to be
%  validated.
%
% Output:
%
%   The results of the validation will be reflected in a sub-field of the
%   photoBleachCorrection field called "validation". If the fit was
%   validated, validation.status will be true, and it will be false
%   otherwise.
% 
% Hunter Elliott
% 12/2009
% 

%% ------ Input ---- %%

movieData = setupMovieData(movieData);

[batchMode,iChannel,doPort,doNorm,alpha] = parseInput(varargin);

% --- Defaults --- %

if isempty(batchMode)
    batchMode = false;
end

if isempty(iChannel) 
    if ~batchMode
        iChannel = selectMovieChannels(movieData,false,'Select channel to validate photobleach correction:');
    else
        error('In batch mode, you must specify the ChannelIndex option!')        
    end
end

if isempty(doPort)
    doPort = true;
end
if isempty(doNorm)
    doNorm = true;
end
if isempty(alpha)
    alpha = .05;
end

%% ------- Init ----- %%

if ~(doPort || doNorm)
    error('At least ONE test must be enabled to validate photobleach correction! Set the PortTest and/or NormTest options to true!')
end

if ~checkMovieProcedure(movieData,'photobleachCorrection',iChannel);
    error('Photobleach correction must be performed before it can be validated! Check photobleach correction!')
end

load([movieData.photobleachCorrection.directory filesep movieData.photobleachCorrection.fileName])


%Check that the needed variables were loaded
if ~exist('resFit','var')
    error('The specified photobleach correction file does not contain the correct variables! Check the photobleach correction!')
end

%% ----- Validation ------ %%

if doPort
    
    maxLag = floor(nnz(~isnan(resFit)) / 4); %Determine the maximum lag based on the t.s. length
    
    hPort = portmanteau(resFit,maxLag,alpha);    
    
else
    hPort = 1;
end

if doNorm
    
    hNorm = kstest(resFit,[],alpha);
   
else
   hNorm = 1; 
end


%% ----- Output ---- %%

if ~hNorm && ~hPort
    movieData.photobleachCorrection.validation.status = 1;
else
    movieData.photobleachCorrection.validation.status = 0;
end

movieData.photobleachCorrection.validation.dateTime = datestr(now);
movieData.photobleachCorrection.validation.iChannel = iChannel;
movieData.photobleachCorrection.validation.portTest = doPort;
movieData.photobleachCorrection.validation.normtest = doNorm;

updateMovieData(movieData);

function [batchMode,iChannel,doPort,doNorm,alpha] = parseInput(argArray)


%Init output
doPort = [];
doNorm = [];
iChannel = [];
batchMode = [];
alpha = [];

if isempty(argArray)
    return
end

nArg = length(argArray);

%Make sure there is an even number of arguments corresponding to
%optionName/value pairs
if mod(nArg,2) ~= 0
    error('Inputs must be as optionName / value pairs!')
end

for i = 1:2:nArg
    
   switch argArray{i}    
       
       case 'BatchMode'
           batchMode = argArray{i+1};
              
       case 'PortTest'
           doPort = argArray{i+1};
           
       case 'NormTest'
           doNorm = argArray{i+1};
           
       case 'ChannelIndex'
           iChannel = argArray{i+1};
           
       case 'Alpha'
           alpha = argArray{i+1};            
           
       otherwise
       
           error(['"' argArray{i} '" is not a valid option name! Please check input!'])
   end
   
end
