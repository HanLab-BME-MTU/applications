function [movieArray,errArray] = batchProcessBiosensorMovies(movieArray,shadeImDirs,actChan,volChan,BTChans,BTcoef,rfpChan)
%BATCHPROCESSBIOSENSORMOVIES sloppy, GUI-free function to processes a series of biosensor movies
%
% All movies must have the same number and ordering of channels. They will
% be processed using the same parameters. This function is horrible, and
% you probably shouldn't use it.
%
% movieArray = batchProcessBiosensorMovies(movieArray,actChan,volChan,shadeImDirs,BTChans,BTcoef)
%
% Hunter Elliott
% 8/2010
%

if nargin < 1 || isempty(movieArray) || ~isa(movieArray(1),'MovieData')
    error('The first input must be an array of MovieData objects!')
end

nMov = numel(movieArray); % Get number of movies to process

if nargin < 2 || isempty(shadeImDirs)
    error('You must input the shade image directories!');
end
if nargin < 3 || isempty(actChan)
    actChan = selectMovieChannels(movieArray(1),0,'Select the activity channel:');
end

if nargin < 4 || isempty(volChan)
    volChan = selectMovieChannels(movieArray(1),0,'Select the volume channel:');
end

if nargin < 5
    BTChans = [];
    disp('No bleedthrough channels specified, skipping BT correction...')
end

if nargin < 6
    rfpChan = [];
end

%% --------- Init ---------- %%

%Common parameters
if ~isempty(BTChans)
    p.ChannelIndex = [actChan volChan rfpChan BTChans(2)];
else
    p.ChannelIndex = [actChan volChan];
end

p.BatchMode = true;

errArray(1:nMov) = struct('hasError',[],'error',[]);


pS = p;      
if ~isempty(BTChans)
    pS.ChannelIndex = [actChan volChan rfpChan BTChans(2)];
else
    pS.ChannelIndex = [actChan volChan rfpChan];
end

pS.ShadeImageDirectories = shadeImDirs;

%% ------------ Processing -------------- %%
%Go through each movie and call the processing function.

%wtBar = waitbar(0,['Please wait, processing movie 1 of ' num2str(nMov)]);
tic;
parfor iMov = 1:nMov
    
    
    try
               
%         pT = p;
%         pT.MaxJump = .25;
%         movieArray(iMov) = thresholdMovie(movieArray(iMov),pT);
%         
        movieArray(iMov) = createMovieBackgroundMasks(movieArray(iMov),p);
                        
        movieArray(iMov) = shadeCorrectMovie(movieArray(iMov),pS);
        
        pBS = p;
        pBS.MaskChannelIndex = p.ChannelIndex;
        movieArray(iMov) = backgroundSubtractMovie(movieArray(iMov),pBS);
        
        if ~isempty(BTChans)
            pBT = p;
            pBT.ChannelIndex = actChan;
            pBT.BleedChannelIndex = BTChans;
            pBT.BleedCoefficients = BTcoef;
            movieArray(iMov) = bleedthroughCorrectMovie(movieArray(iMov),pBT);            
        end
% 
%         pR = p;
%         pR.ChannelIndex = [actChan volChan];
%         movieArray(iMov) = ratioMovie(movieArray(iMov),pR);
%         
%         pPB = p;        
%         pPB.ChannelIndex = actChan;        
%         movieArray(iMov) = photobleachCorrectMovieRatios(movieArray(iMov),pPB);
%         
    
    catch errMess

        errArray(iMov).hasError = true;
        errArray(iMov).error = errMess;
        
        disp(['======= Error processing movie ' num2str(iMov) ' : ' errMess.message '======']);            
    
    end
    
    disp(['**** FINISHED MOVIE ' num2str(iMov) ' OF ' num2str(nMov) '*****'])
    
    %waitbar(iMov/nMov,wtBar,['Please wait, processing movie ' num2str(iMov) ' of ' num2str(nMov)]);
      
    
end

toc

% if ishandle(wtBar) 
%     close(wtBar)
% end

