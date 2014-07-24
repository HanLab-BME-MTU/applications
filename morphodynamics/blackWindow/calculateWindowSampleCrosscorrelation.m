function movieData = calculateWindowSampleCrosscorrelation(movieData,paramsIn)
%CALCULATEWINDOWSAMPLEAUTOCORRELATION calcualates and makes figures of the crosscorrelation of the window samples 
% 
% movieData = calculateWindowSampleAutocorrelation(movieData)
% movieData = calculateWindowSampleAutocorrelation(movieData)
% 
% This function calculates the crosscorrelation between the window samples
% and protrusion samples for each window of the input movie. Window
% sampling must already have been run using sampleMovieWindows.m. The
% averaged crosscorrelation along the cell edge, into the cell, etc are
% also calcualted and figures created for each. Additionally, the average
% cross-correlation of a window with its neighbors is calculated and
% returned.
% 
% Input:
% 
%   movieData - The MovieData object describing the movie to analyze.
% 
%   paramsIn - A structure containing the parameters to use when
%   calcualating crosscorrelatioins. The parameters should be stored as
%   fields in the structure, with the field names and possible values as
%   described below:
% 
%   Possible Parameter Structure Field Names:
%       ('FieldName' -> possible values)
%
%       ('MaxLag'->positive integer scalar) The maximum lag, in frames, to
%       calculate the crosscorrelation to. Optional. Default is 1/4 of the
%       total number of frames.
%
%       ('UseBands' -> Positive integer scalar/Vector) Specifies which
%       bands of windows to use in calculating combined cross-correlation.
%       Optional. Default is to use all bands.
%
%       ('OutputDirectory' -> character string) Optional. A character
%       string specifying the directory to save the results and figures to.
%       Samples for different channels will be saved as files in this
%       directory.
%       If not input, the results will be saved to the same directory as the
%       movieData, in a sub-directory called "window_sample_crosscorrelation"
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

%% ------------------- Params ------------------ %%

%Number of lags must be <= N/4, N must be >50
%Ref: Time Series Analysis, Forecast and Control. Jenkins, G. Box,G
nPmin = 50;

%% ----------------- Input -------------- %%



if nargin < 1 || ~isa(movieData,'MovieData')   
    error('The first input must be a valid MovieData object!');        
end

if nargin < 2 || isempty(paramsIn)
    paramsIn = [];  
end

%Check if crosscorrelation has been run before 
iProc = movieData.getProcessIndex('SampleCrosscorrelationProcess',1,false);
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(SampleCrosscorrelationProcess(movieData,movieData.outputDirectory_));
end

%Parse input, store in parameter structure
p = parseProcessParams(movieData.processes_{iProc},paramsIn);

%Make sure the movie has been sampled
iSampProc = movieData.getProcessIndex('WindowSamplingProcess',1,false);
if isempty(iSampProc) || ~movieData.processes_{iSampProc}.checkChannelOutput(p.ChannelIndex)
    error('The movie must have valid window sampling for all selected channels! Please run sampleMovieWindows.m first!')   
end

%Get the windowing process for later use
%iWinProc = movieData.getProcessIndex('WindowingProcess',1,0);

%Make sure the protrusion calculation has been run
iProtProc = movieData.getProcessIndex('ProtrusionSamplingProcess',1,~p.BatchMode);
if isempty(iProtProc) || ~movieData.processes_{iProtProc}.checkChannelOutput
    disp('Cannot calculate protrusion-cross correlation: No protrusion calculation found!')
else    
    protSamples = movieData.processes_{iProtProc}.loadChannelOutput;    
end

if movieData.nFrames_ < nPmin
    error('Too few frames to calculate cross-correlation!!! Need at least 50 frames!')
end



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
    tData = -p.MaxLag*movieData.timeInterval_:movieData.timeInterval_:p.MaxLag*movieData.timeInterval_;
    tLabel = 'Time Lag, Seconds';
else
    tData = -p.MaxLag:p.MaxLag;
    tLabel = 'Time Lag, Frames';
end
   
%Make sure the max lag isn't too large for the data.
p.MaxLag = min(p.MaxLag,floor(movieData.nFrames_/4));


%Double-exponential function for photobleach correction fit
fitFun = @(b,x)(b(1) .* exp(b(2) .* x))+(b(3) .* exp(b(4) .* x));
bInit = [1 0 1 0]; %Initial guess for fit parameters.
fitOptions = statset('Robust','on','MaxIter',5e3,'Display','off');
timePoints = 0:(movieData.nFrames_-1);

if isempty(p.UseBands)
    bandsSelected = false;
else
    bandsSelected = true;    
    nBandUse = numel(p.UseBands);
end

%% ----------------- Crosscorr Calc ------------------ %%

for iChan = 1:nChan
    
    %Load the activity samples for this channel
    actSamples = movieData.processes_{iSampProc}.loadChannelOutput(p.ChannelIndex(iChan));        
    
    [nStripMax,nBandMax,~] = size(actSamples.avg);    
    
    if ~bandsSelected
        p.UseBands = 1:nBandMax;
        nBandUse = numel(p.UseBands);
    end
    
    nObsAct = zeros(nStripMax,nBandMax);
    nObsProt = zeros(nStripMax,1);
    nObsComb = zeros(nStripMax,nBandMax);
        
    ccProtActPerWin = nan(nStripMax,nBandMax,2*p.MaxLag+1,2);
        
    isActStationary = false(nStripMax,nBandMax);
    isProtStationary = false(nStripMax,1);
    
    %Photobleach correction.
    %Fit function to ratio timeseries  TEMP - this needs to be optional in
    %case correction has already been performed!!!!!!!!!!!!!!!!!!!!!!!!

    fitData = squeeze(nanmean(nanmean(actSamples.avg,1),2));
    [bFit,resFit,jacFit,covFit,mseFit] = nlinfit(timePoints(:),fitData,fitFun,bInit,fitOptions);
    %Get confidence intervals of fit and fit values
    [fitValues,deltaFit] = nlpredci(fitFun,timePoints(:),bFit,resFit,'covar',covFit,'mse',mseFit);
    
    
    %Check the fit jacobian
    [dummy,R] = qr(jacFit,0); %#ok<ASGLU>
    if condest(R) > 1/(eps(class(bFit)))^(1/2)        
        warning('WARNING: The photobleach correction fit is not very good. This may intruduce artifacts in the resulting cross-correlation!!') %#ok<WNTAG>
    end
       
    if p.BatchMode
        pbCorrFig = figure('Visible','off');
    else
        pbCorrFig = figure;
    end    
    hold on
    title('Photobleach correction fit')
    xlabel('Frame Number')
    plot(timePoints,fitData)
    plot(timePoints,fitValues,'r')
    plot(timePoints,fitValues+deltaFit,'--r')
    legend('Average Activity','Fit','Fit 95% C.I.')
    plot(timePoints,fitValues-deltaFit,'--r')

    hgsave(pbCorrFig,[p.OutputDirectory filesep 'photobleach correction fit.fig']);
    
    for k = 1:nStripMax 
        
        nObsProt(k) = nnz(~isnan(protSamples.avgNormal(k,:)));
        if nObsProt(k) > nPmin
            isProtStationary(k) = adftest(squeeze(protSamples.avgNormal(k,:)));
        end
        for l = 1:nBandMax
                        
            %Get number of non-Nan points
            nObsAct(k,l) = nnz(~isnan(actSamples.avg(k,l,:)));
            
            %Apply the photobleach correction
            actSamples.avg(k,l,:) = squeeze(actSamples.avg(k,l,:)) ./ fitValues;
            
            if nObsAct(k,l) > nPmin                                
                
                %Test the time-series for stationarity using Dickey-Fuller.
                %Uses default alpha of .05;
                isActStationary(k,l) = adftest(squeeze(actSamples.avg(k,l,:)) - nanmean(squeeze(actSamples.avg(k,l,:))));           

                %Calculate the cross-corr for this window, if there are
                %enough observations and at least one of the time series is
                %stationary.
                if nObsAct(k,l) >= nPmin && nObsProt(k) >= nPmin && isProtStationary(k) && isActStationary(k,l)
                    ccProtActPerWin(k,l,:,:) = crossCorr(squeeze(actSamples.avg(k,l,:)),protSamples.avgNormal(k,:)',p.MaxLag);
                    nObsComb(k,l) = nnz(squeeze(~isnan(actSamples.avg(k,l,:))) & ~isnan(protSamples.avgNormal(k,:))');
                    if any(abs(ccProtActPerWin(k,l,:,1))>1) %TEMP? - workaround for the fact that even with 50 points, we can have CC values outside the range {-1, 1}, presumably because the STD estimate can be horrible depending on the NaN pattern
                        ccProtActPerWin(k,l,:,:) = NaN;
                    end
                end                                   
            end
        end                
            
    end                                    
    
    if p.BatchMode
        ccFig = figure('Visible','off');
    else
        ccFig = figure;
    end
    
    %Bootstrap the cross-correlation for the selected bands    
    
    %Extract and reshape the individual correlations to be combined
    allCC = reshape(ccProtActPerWin(:,p.UseBands,:,1),nStripMax*nBandUse,2*p.MaxLag+1)';    
    allCB = 1.96 ./ sqrt(reshape(nObsComb(:,p.UseBands),nStripMax*nBandUse,1)');
    hasCorr = sum(~isnan(allCC),1) > 0;
    allCC = allCC(:,hasCorr);
    allCB = allCB(hasCorr);
    
    [combMeanCC,combBootCI] = correlationBootstrap(allCC,allCB);        

    plot(tData,combMeanCC)
    hold on    
    plot(tData,combBootCI(1,:),'--')
    legend('Mean Correlation','Bootstrapped 95% CI')
    plot(tData,combBootCI(2,:),'--')
    plot([0 0],ylim,'--k')
    plot(xlim,[0 0],'--k')
    xlabel(tLabel); %TEMP - get time interval and scale if available!!!!
    ylabel('Cross Correlation');
    %TEMP - Change this so that if all bands are used, it says "all
    %samples" and if not, it says the bands used
    title({['Selected Bands, All Strips Combined, Temporal Cross-Correlation of Protrusion and Activity, Channel ' num2str(p.ChannelIndex(iChan))],'Positive lags mean protrusion follows activity'});
    
    hgsave(ccFig,[p.OutputDirectory filesep 'combined temporal crosscorrelation between protrusion and activity channel ' num2str(p.ChannelIndex(iChan))]);            
            
    
    %Show the individual contributions to the combined cross-corr   
    if p.BatchMode
        ccFigBand = figure('Visible','off');
        ccFigBandMean = figure('Visible','off');
    else
        ccFigBand = figure;
        ccFigBandMean = figure;
    end         
    for j = 1:nBandUse
        figure(ccFigBand)
        subplot(1,nBandUse,j);
        hold on;
        imagesc(tData,1:nStripMax,squeeze(ccProtActPerWin(:,p.UseBands(j),:,1)))
        xlabel(tLabel),ylabel('Window #, Along cell edge'),colorbar
        plot([0 0],ylim,'--w')
        title(['Per-Strip Cross-correlation, band ' num2str(p.UseBands(j))])
        
        figure(ccFigBandMean)
        subplot(1,nBandUse,j);
        hold on;
        bandCC = squeeze(ccProtActPerWin(:,p.UseBands(j),:,1))';        
        bandCB = 1.96 ./ sqrt(nObsComb(:,p.UseBands(j)))';
        hasCorr = sum(~isnan(bandCC),1) > 0;
        bandCC = bandCC(:,hasCorr);
        bandCB = bandCB(hasCorr);
        [bandMeanCC(j,:),bandBootCI(j,:,:)] = correlationBootstrap(bandCC,bandCB);                          
        plot(tData,bandMeanCC(j,:))
        plot(tData,squeeze(bandBootCI(j,1,:)),'--')
        legend('Mean CrossCorr','Bootstrappedn 95% CI')
        plot(tData,squeeze(bandBootCI(j,2,:)),'--')
        plot([0 0 ],ylim,'--k');
        plot(xlim,[ 0 0 ],'--k');
        title(['All Strip mean correlation, band ' num2str(p.UseBands(j))])
        xlabel(tLabel);
        ylabel('Cross Correlation')
        
    end
    
    hgsave(ccFigBand,[p.OutputDirectory filesep 'per window temporal crosscorrelation between protrusion and activity channel ' num2str(p.ChannelIndex(iChan))]);            
    hgsave(ccFigBandMean,[p.OutputDirectory filesep 'per band mean temporal crosscorrelation between protrusion and activity channel ' num2str(p.ChannelIndex(iChan))]);                                
    
    save([p.OutputDirectory filesep 'temporal crosscorrelation channel ' num2str(p.ChannelIndex(iChan))],'ccProtActPerWin','nObsComb','nObsAct','nObsProt','combMeanCC','combBootCI','bandMeanCC','bandBootCI');        
            
end


%Set output directories in movieData
%TEMP -DOOOOOO THHHIIISSS!!!!




