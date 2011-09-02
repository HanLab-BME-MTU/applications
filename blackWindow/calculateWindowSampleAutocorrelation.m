function movieData = calculateWindowSampleAutocorrelation(movieData,paramsIn)
%CALCULATEWINDOWSAMPLEAUTOCORRELATION calcualates and makes figures of the autocorrelation of the window samples 
% 
% movieData = calculateWindowSampleAutocorrelation(movieData)
% movieData = calculateWindowSampleAutocorrelation(movieData)
% 
% This function calculates the autocorrelation for the samples for each
% window of the input movie. Window sampling must already have been run
% using sampleMovieWindows.m. The averaged autocorrelation along the cell
% edge, into the cell, etc are also calcualted and figures created for
% each.
% 
% Input:
% 
%   movieData - The MovieData object describing the movie to analyze.
% 
%   paramsIn - A structure containing the parameters to use when
%   calcualating autocorrelatioins. The parameters should be stored as
%   fields in the structure, with the field names and possible values as
%   described below:
% 
%   Possible Parameter Structure Field Names:
%       ('FieldName' -> possible values)
%
%       ('MaxLag'->positive integer scalar) The maximum lag, in frames, to
%       calculate the autocorrelation to. Optional. Default is 1/4 of the
%       total number of frames.
%
%       ('DetrendMethod'->'linear','robustLinear','difference','none') This
%       specifies what kind of de-trending to apply to the samples prior to
%       calculating autocorrelation:
%           'linear' - Default. Linear least-squares fit is removed from data.
%           'robustLinear' - Linear least-median-squares fit is removed
%              from data.
%           'difference' - The first difference is taken.
%           'none' - What do you think this one does?
%
%       ('OutputDirectory' -> character string) Optional. A character
%       string specifying the directory to save the results and figures to.
%       Samples for different channels will be saved as files in this
%       directory.
%       If not input, the results will be saved to the same directory as the
%       movieData, in a sub-directory called "window_sample_autocorrelation"
%
%       ('ChannelIndex' -> Positive integer scalar or vector) Optional. The
%       integer index of the channels to analyze samples from. If not input,
%       all channels which have been sampled will be used.
%
%       ('BatchMode' -> True/False)
%       If true, graphical output, figure display and user interaction is
%       supressed (i.e. progress bars, dialog and question boxes etc.) Note
%       that if this is set to true, the .fig versions of the figures will
%       be saved as invisible and you will have to set the 'Visible'
%       property to 'On' after opening them to see them.
%
%
% Output:
%
%   movieData - The updated MovieData object with the location of the
%   results logged in the processes_ array. Additionally, the results and
%   figures will be saved to the specified OutputDirectory.
%
% Hunter Elliott
% 8/2011
%

%% -------------------------- Input --------------------------------%%

if nargin < 1 || ~isa(movieData,'MovieData')   
    error('The first input must be a valid MovieData object!');        
end

if nargin < 2 || isempty(paramsIn)
    paramsIn = [];  
end

%Check if autocorrelation has been run before 
iProc = movieData.getProcessIndex('SampleAutocorrelationProcess',1,false);
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(SampleAutocorrelationProcess(movieData,movieData.outputDirectory_));
end

%Parse input, store in parameter structure
p = parseProcessParams(movieData.processes_{iProc},paramsIn);

%Make sure the movie has been sampled
iSampProc = movieData.getProcessIndex('WindowSamplingProcess',1,false);
if isempty(iSampProc) || ~movieData.processes_{iSampProc}.checkChannelOutput(p.ChannelIndex)
    error('The movie must have valid window sampling for all selected channels! Please run sampleMovieWindows.m first!')   
end

%Get the windowing process for later use
iWinProc = movieData.getProcessIndex('WindowingProcess',1,0);

%% ---------------------- Init ---------------------------- %%

%Set up the output directory
mkClrDir(p.OutputDirectory);

nChan = numel(p.ChannelIndex);

%Print options for creating .eps files.
pOpt = {'-r300',...% dpi = 300
        '-depsc2',...% use eps format
                };

%Get time data if available
if ~isempty(movieData.timeInterval_)
    tData = 0:movieData.timeInterval_:p.MaxLag*movieData.timeInterval_;
    tLabel = 'Time Lag, Seconds';
else
    tData = 0:p.MaxLag;
    tLabel = 'Time Lag, Frames';
end
            

%% ------------------ Autocorr Calc ----------------------- %%

disp('Calculating window sample autocorrelation....')

for iChan = 1:nChan    

    samples = movieData.processes_{iSampProc}.loadChannelOutput(p.ChannelIndex(iChan));
    
    [nStripMax,nBandMax,nFrames] = size(samples.avg);
    
    actTimeSeries(1:nStripMax,1:nBandMax) = struct('observations',[]);
    nObs = nan(nStripMax,nBandMax);
    acPerWin = nan(nStripMax,nBandMax,p.MaxLag+1,2);
    acPerBand = nan(p.MaxLag+1,2,nBandMax);
    acPerStrip = nan(p.MaxLag+1,2,nStripMax);            

    for l = 1:nBandMax

        for k = 1:nStripMax
            %Extract the current time-series
            actTimeSeries(k,l).observations = squeeze(samples.avg(k,l,:));

            %Count the observations, leave the empty series as placeholders
            nObs(k,l) = nnz(~isnan(actTimeSeries(k,l).observations));                

            %De-trend the current time-series.
            if nObs(k,l) >= 3 %We need at least 3 points for de-trending to make sense                    
                switch p.DetrendMethod
                    
                    case 'linear'
                        actTimeSeries(k,l).observations = removeLinearTrend(actTimeSeries(k,l).observations,[],0);                        
                    case 'linearRobust'
                        actTimeSeries(k,l).observations = removeLinearTrend(actTimeSeries(k,l).observations,[],1);
                    case 'difference'
                        actTimeSeries(k,l).observations = diff(actTimeSeries(k,l).observations);
                    case 'none'
                        
                    otherwise
                        error(['"' p.DetrendMethod '" is not a supported detrending method! Check the DetrendMethod option!'])
                end                
                %If its long enough, get the individual autocorrelation
                if (nObs(k,l) - p.MaxLag) >= (3*p.MaxLag)
                    acPerWin(k,l,:,:) = autoCorr(actTimeSeries(k,l),p.MaxLag,-1);
                end
            end

        end
        %Get the per-band average autocorrelation
        if (sum(squeeze(nObs(nObs(:,l)>p.MaxLag,l))) - p.MaxLag) >= (3*p.MaxLag) %Make sure we have enough points
            acPerBand(:,:,l) = autoCorr(actTimeSeries(nObs(:,l)>=3,l),p.MaxLag,-1);%we need at least 3 points per traj for trend removal
        end

    end                
    %Get the per-strip average autocorrelation
    for k = 1:nStripMax
        if (sum(squeeze(nObs(k,nObs(k,:)>p.MaxLag))) - p.MaxLag) >= (3*p.MaxLag) %Make sure we have enough points by K's criteria
            acPerStrip(:,:,k) = autoCorr(actTimeSeries(k,nObs(k,:)>=3),p.MaxLag,-1);
        end
    end            
    
    nPtsTot = sum(nObs(nObs(:)>=3));
    %Get the all-sample combined autocorrelation
    [acActTimeAll,~] = autoCorr(actTimeSeries(nObs(:)>=3),p.MaxLag,-1);

    if p.BatchMode
        acFig = figure('Visible','off');
    else
        acFig = figure;
    end

    plot(tData,acActTimeAll(:,1),'r');
    hold on
    plot(xlim,ones(1,2)*1.96/sqrt(nPtsTot),'--k');
    plot(xlim,-ones(1,2)*1.96/sqrt(nPtsTot),'--k');
    plotTransparent(tData,acActTimeAll(:,1),acActTimeAll(:,2),'r',.2,0);
    xlabel(tLabel);
    ylabel('Autocorrelation');
    title(['All Samples Combined, Temporal Autocorrelation, Channel ' num2str(p.ChannelIndex(iChan))]);
    
    currFileName = [p.OutputDirectory filesep 'combined temporal autocorrelation channel ' num2str(p.ChannelIndex(iChan))];
    hgsave(acFig,currFileName);            
    print(pOpt{:},currFileName);
    
    
    
    if p.BatchMode
        acFigBand = figure('Visible','off');
    else
        acFigBand = figure;
    end       
    
    surf(squeeze(acPerBand(:,1,:)))
    view([-125 20])
    xlabel('Into Cell (band #)')
    ylabel('Time Lag (frames)')
    zlabel('Autocorrelation')
    
    title(['Per Band, Temporal Autocorrelation, Channel ' num2str(p.ChannelIndex(iChan))]);
    
    currFileName = [p.OutputDirectory filesep 'per band temporal autocorrelation channel ' num2str(p.ChannelIndex(iChan))];
    hgsave(acFigBand,currFileName);
    print(pOpt{:},currFileName);
    
    
    if p.BatchMode
        acFigStrip = figure('Visible','off');
    else
        acFigStrip = figure;
    end            
    
    surf(squeeze(acPerStrip(:,1,:)))
    view([-125 20])
    xlabel('Along Cell Edge (strip #)')
    ylabel('Time Lag (frames)')
    zlabel('Autocorrelation')
    
    title(['Per Strip, Temporal Autocorrelation, Channel ' num2str(p.ChannelIndex(iChan))]);
    
    currFileName = [p.OutputDirectory filesep 'per strip temporal autocorrelation channel ' num2str(p.ChannelIndex(iChan))];
    hgsave(acFigStrip,currFileName);        
    print(pOpt{:},currFileName);
    
    if p.BatchMode
        acFigWin = figure('Visible','off');
    else
        acFigWin = figure;
    end    
    acSum = sqrt(nansum(acPerWin(:,:,:,1) .^2,3));
    nCol = 1024;
    winColors = colormap(jet(nCol));
    hold on
    %maxAc = mean(acSum(acSum(:)>0))+std(acSum(acSum(:)>0));
    maxAc = max(acSum(:));
    minAc = 1;    
        
    windows = movieData.processes_{iWinProc}.loadChannelOutput(1);    
    nStrip = numel(windows);
    for k = 1:nStrip
        nBandCur = numel(windows{k});
        for l = 1:nBandCur

            if acSum(k,l) > 0
                iCol = round(interp1(linspace(minAc,maxAc,nCol),1:nCol,acSum(k,l)));%This is probably slow - do it manually?                    
                plotString = {winColors(iCol,:)};
            else
                plotString = {'k','FaceAlpha',0};
            end
            if ~isempty(windows{k}{l})
                plotWindows(windows{k}{l},plotString)
            end                
        end
    end    
            
    axis off
    caxis([minAc maxAc])
    colorbar
    title('Root-Sum-Squared Autocorrelation at all Lags, Per window')
    
    currFileName = [p.OutputDirectory filesep 'per window temporal autocorrelation channel ' num2str(p.ChannelIndex(iChan))];
    hgsave(acFigWin,currFileName);
    print(pOpt{:},currFileName);
                
    save([p.OutputDirectory filesep 'temporal autocorrelation channel ' num2str(p.ChannelIndex(iChan))],'acActTimeAll','acPerBand','acPerStrip','acPerWin','nObs');   
    
end

    



