function makeWindowMovie(movieData,varargin)
%MAKEWINDOWMOVIE makes an avi or .mov movie with windows overlain on images from the input movie
%
% movieData = makeMaskMovie(movieData)
% 
% movieData = makeMaskMovie(movieData,'OptionName',OptionValue)
%
% Makes a movie overlaying the sampling windows on the images from each
% channel so that the user can verify which areas of the image are being
% sampled.
%
%
% Input:
%
%   movieData - The MovieData object describing the movie to make a window
%   movie for, as created with setupMovieDataGUI.m and windowed with
%   getMovieWindows.m
%
%   Possible Option Names:
%
%   ('OptionName'->Possible values)
%
%       ('ChannelIndex'-> positive integer scalar or vector of length <=3)
%       These numbers correspond to the indices of the channel(s) to overlay
%       windows on. 
%
%       ('FigureHandle' - Positive integer) figure handle to plot the
%       images / windows on.
%       Optional. Default is to create a new figure.
%
%       ('FileName'-> Character array) String specifying file name to save
%       movie as.
%       Optional. Default is "windowMovie"
%
%       ('MakeAvi' - Logical scalar) If true, the movie will be saved as .avi.
%       Optional. Default is false.
%
%       ('MakeMov' - Logical scalar) If true, movie will be saved as .mov.
%       Optional. Default is true.
%
% Output:
%
%   movieData - The movieData with the movie's name and location logged in
%   it.
%
%   The resulting movie will be saved in the analysis directory specified in
%   the movieData.
% 
% Hunter Elliott
% 8/2010
%       Optional. If not input, user is asked.
%

%% ------- Input ------------ %%


% Parse input from variable input arguments
[iChannels,figHan,mvName,makeAvi,makeMov] = parseInput(varargin);

%Make sure that either makAvi or makeMov are true!!
if ~(makeMov || makeAvi)
    error('Either makeMov or makeAvi MUST be true!! (otherwise no movie will be saved!!)')
end

%Check/set up the movieData
if ~isa(movieData,'MovieData');
    error('The first input must be a valid MovieData object!');
end

if isempty(figHan)    
    figHan = figure;    
elseif ishandle(figHan)
    figure(figHan)    
else
    error('Input figure handle is not a valid graphics handle!!')
end

%Make sure the movie has been segmented
iWinProc = movieData.getProcessIndex('WindowingProcess',1,true);

if isempty(iWinProc)
    error('The movie must have been windowed before you can make a window overlay movie!')    
end

if isempty(iChannels)
    iChannels = selectMovieChannels(movieData,true);
end

nChanMov = length(iChannels);

if nChanMov == 0 || nChanMov > 3
    error('The numer of channels selected must be between 1 and 3 !!')
end


%% ------ Init ------ %%


nImages = movieData.nFrames_;
figure(figHan);
axHandle = gca;



%% ------ Movie making -----%%

%Go through the frames/channels and make the movie
for j = 1:nImages    
    
    %Show the image+windows
    imageViewer(movieData,'ChannelIndex',iChannels,...
                           'Frame',j,...                           
                           'Overlay','Windows','AxesHandle',axHandle);
            
    if makeAvi
        windowMovie(j) = getframe(figHan); %#ok<AGROW> Pre-allocation not needed with movies. I'm serious.
    end    
           
    if makeMov        
        if j == 1
            MakeQTMovie('start',[movieData.outputDirectory_ filesep mvName '.mov'])
            MakeQTMovie('quality',.85)
        end
        MakeQTMovie('addfigure')
        
    end
    cla(axHandle)
end

%% ----- Finalization ----- %%

if makeMov
    MakeQTMovie('finish')        
end

if makeAvi
    %Make the movie, save it to the movie directory
    if isunix %Compression is not supported under unix!
        movie2avi(windowMovie,[movieData.outputDirectory_  filesep mvName '.avi']);        
    else
        movie2avi(windowMovie,[movieData.outputDirectory_  filesep mvName '.avi'],'compression','Cinepak')    
    end    
end

if ishandle(figHan)%Make sure the user hasn't closed it already.
    close(figHan);
end

function [iChannels,figHan,mvName,makeAvi,makeMov] = parseInput(argArray)
%Sub-function for parsing variable input arguments



%-----Defaults-----%
iChannels = [];
figHan = [];
mvName = 'windowMovie';
makeAvi = false;
makeMov = true;

if isempty(argArray)
    return
end

nArg = length(argArray);

%Make sure there is an even number of arguments corresponding to
%optionName/value pairs
if mod(nArg,2) ~= 0
    error('Inputs must be as optionName/ value pairs!')
end

for i = 1:2:nArg
    
   switch argArray{i}                     
              
       case 'ChannelIndex'           
           iChannels = argArray{i+1};
           
       case 'FigureHandle'           
           figHan = argArray{i+1};
              
       case 'FileName'
           mvName = argArray{i+1};
           
       case 'MakeAvi'
           makeAvi = argArray{i+1};
           
       case 'MakeMov'
           makeMov = argArray{i+1};
           
           
       otherwise
       
           error(['"' argArray{i} '" is not a valid option name! Please check input!'])
   end
               
      
   
end
