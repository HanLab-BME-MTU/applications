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
%       ('MaxBands' ->Positive integer scalar) The maximum number of bands
%       of windows, starting at the cell edge, to include in the combined
%       cross-correlation. Because the correlations with edge velocity
%       often fall off further from the edge, specifying a smaller number
%       here may give a stronger combined cross-correlation. Optional.
%       Default is Infinity - that is, all bands are included.
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
            
   
%% ----------------- Crosscorr Calc ------------------ %%

for iChan = 1:nChan
    
    %Load the activity samples for this channel
    actSamples = movieData.processes_{iSampProc}.loadChannelOutput(p.ChannelIndex(iChan));        
    
    [nStripMax,nBandMax,~] = size(actSamples.avg);
    
    nObs = zeros(nStripMax,nBandMax);
    
    actTimeSeries(1:nStripMax,1:nBandMax) = struct('observations',[]);
    protTimeSeries(1:nStripMax) = struct('observations',[]);
    ccProtActPerWin = nan(nStripMax,nBandMax,2*p.MaxLag+1,2);
    ccProtActPerStrip = nan(nStripMax,2*p.MaxLag+1,2);
    ccProtActPerBand = nan(nBandMax,2*p.MaxLag+1,2);
    
    %Get the current maximum bands
    currMaxBands = min(p.MaxBands,nBandMax);
    
    for k = 1:nStripMax 
        
        for l = 1:nBandMax
            
            %Get the protrusion for this strip. We copy this into every
            %band to make the combined cross-correlation easier later.
            protTimeSeries(k,l).observations = squeeze(protSamples.avgNormal(k,:))';
        

            %Extract the current activity time-series
            actTimeSeries(k,l).observations = squeeze(actSamples.avg(k,l,:));
            
            %Leave the empty series as placeholders, and get number of
            %non-Nan points
            nObs(k,l) = nnz(~isnan(actTimeSeries(k,l).observations));
            
            %Calculate the cross-corr for this window
            if (nObs(k,l) - p.MaxLag) >= 3*p.MaxLag                                 
                ccProtActPerWin(k,l,:,:) = crossCorr(actTimeSeries(k,l),protTimeSeries(k,l),p.MaxLag);                                                                          
            end                                   
            
        end
        
        if (sum(squeeze(nObs(k,nObs(k,1:currMaxBands)>p.MaxLag))) - p.MaxLag) >= (3*p.MaxLag) %Make sure we have enough points
            ccProtActPerStrip(k,:,:) = crossCorr(actTimeSeries(k,1:currMaxBands),protTimeSeries(k,1:currMaxBands),p.MaxLag);
        end
            
    end                
        
    for l = 1:nBandMax
        if (sum(squeeze(nObs(nObs(:,l)>p.MaxLag))) - p.MaxLag) >= (3*p.MaxLag) %Make sure we have enough points by K's criteria
            ccProtActPerBand(l,:,:) = crossCorr(actTimeSeries(nObs(:,l)>=3,l),protTimeSeries(nObs(:,l)>=3,l),p.MaxLag);
        end
    end        

    %Go through again and do local activity cross-corr. Better way to do this than looping through again???
    %Lazy right now....
    ccActLocal = nan(nStripMax,nBandMax,3,3,2*p.MaxLag+1,2);
    for k = 1:nStripMax
        
        for l = 1:nBandMax
            
            if (nObs(k,l) - p.MaxLag) >= 3*p.MaxLag      
                
                %Throw the auto-corr in for good measure
                ccActLocal(k,l,2,2,:,:) = crossCorr(actTimeSeries(k,l),actTimeSeries(k,l),p.MaxLag);                
                
                %Do cross-correlatioin with neighborhing time-series
                %TEMP: If whole-cell (closed contours) this should "wrap
                %around" - HLE
                %Also, there's probably a better way to do
                %this, without all the if statements...????
                if k > 1 && (nObs(k-1,l) - p.MaxLag) >= 3*p.MaxLag
                    ccActLocal(k,l,2,1,:,:) = crossCorr(actTimeSeries(k,l),actTimeSeries(k-1,l),p.MaxLag);                           
                    if l > 1 && (nObs(k-1,l-1) - p.MaxLag) >= 3*p.MaxLag                
                        ccActLocal(k,l,1,1,:,:) = crossCorr(actTimeSeries(k,l),actTimeSeries(k-1,l-1),p.MaxLag);                                                   
                    end
                    if l < nBandMax && (nObs(k-1,l+1) - p.MaxLag) >= 3*p.MaxLag
                        ccActLocal(k,l,3,1,:,:) = crossCorr(actTimeSeries(k,l),actTimeSeries(k-1,l+1),p.MaxLag);                        
                    end                    
                end
                
                if l > 1  && (nObs(k,l-1) - p.MaxLag) >= 3*p.MaxLag                                    
                    ccActLocal(k,l,1,2,:,:) = crossCorr(actTimeSeries(k,l),actTimeSeries(k,l-1),p.MaxLag);                    
                    if k < nStripMax && (nObs(k+1,l-1) - p.MaxLag) >= 3*p.MaxLag                        
                        ccActLocal(k,l,1,3,:,:) = crossCorr(actTimeSeries(k,l),actTimeSeries(k+1,l-1),p.MaxLag);                    
                    end                    
                end
                if l < nBandMax && (nObs(k,l+1) - p.MaxLag) >= 3*p.MaxLag                        
                    ccActLocal(k,l,3,2,:,:) = crossCorr(actTimeSeries(k,l),actTimeSeries(k,l+1),p.MaxLag);                    
                    if k < nStripMax && (nObs(k+1,l+1) - p.MaxLag) >= 3*p.MaxLag                        
                        ccActLocal(k,l,3,3,:,:) = crossCorr(actTimeSeries(k,l),actTimeSeries(k+1,l+1),p.MaxLag);
                    end
                end
                if k < nStripMax && (nObs(k+1,l) - p.MaxLag) >= 3*p.MaxLag                        
                    ccActLocal(k,l,2,3,:,:) = crossCorr(actTimeSeries(k,l),actTimeSeries(k+1,l),p.MaxLag);
                end
                
            end            
        end
    end
    
    
    ccProtActAll = crossCorr(actTimeSeries(:,1:currMaxBands),protTimeSeries,p.MaxLag);       
    
    if p.BatchMode
        ccFig = figure('Visible','off');
    else
        ccFig = figure;
    end

    plot(-p.MaxLag:p.MaxLag,ccProtActAll(:,1),'r');
    hold on    
    plotTransparent(-p.MaxLag:p.MaxLag,ccProtActAll(:,1),ccProtActAll(:,2),'r',.2,0);
    plot([0 0],ylim,'--k')
    xlabel('Time Lag, frames'); %TEMP - get time interval and scale if available!!!!
    ylabel('cross-correlation');
    %TEMP - Change this so that if all bands are used, it says "all
    %samples" and if not, it says the bands used
    title({['Selected Bands, All Strips Combined, Temporal Cross-Correlation of Protrusion and Activity, Channel ' num2str(p.ChannelIndex(iChan))],'Positive lags mean protrusion follows activity'});
    
    hgsave(ccFig,[p.OutputDirectory filesep 'combined temporal crosscorrelation between protrusion and activity channel ' num2str(p.ChannelIndex(iChan))]);            
            
    if p.BatchMode
        ccFigBand = figure('Visible','off');
    else
        ccFigBand = figure;
    end        
    
    imagesc(-p.MaxLag:p.MaxLag,1:nBandMax,squeeze(ccProtActPerBand(:,:,1)))    
    ylabel('Into Cell (band #)')
    xlabel('Time Lag (frames)')
    colorbar
    axis image
    
    title(['Per Band, Temporal Cross-Correaltion of Protrusion and Activity, Channel ' num2str(p.ChannelIndex(iChan))]);
    
    hgsave(ccFigBand,[p.OutputDirectory filesep 'per band temporal crosscorrelation between protrusion and acitivity channel ' num2str(p.ChannelIndex(iChan))]);
    
    if p.BatchMode
        ccFigStrip = figure('Visible','off');
    else
        ccFigStrip = figure;
    end        
    
    imagesc(-p.MaxLag:p.MaxLag,1:nStripMax,squeeze(ccProtActPerStrip(:,:,1)))    
    axis image
    ylabel('Along Cell Edgel (band #)')
    xlabel('Time Lag (frames)')
    colorbar
    
    title(['Per Strip, Selected Bands, Temporal Cross-Correaltion of Protrusion and Activity, Channel ' num2str(p.ChannelIndex(iChan))]);
    
    hgsave(ccFigStrip,[p.OutputDirectory filesep 'per strip temporal crosscorrelation between protrusion and acitivity channel ' num2str(p.ChannelIndex(iChan))]);

    if p.BatchMode
        ccFigLocal = figure('Visible','off');
    else
        ccFigLocal = figure;
    end        
    title(['Cross-Correlation of activity with neighboring windows, Channel ' num2str(p.ChannelIndex(iChan))]);
    hold on
    titlesArray = {'Prev Strip, Towards Edge','Same Strip, Towards Edge','Next Strip, Towards Edge';
                   'Prev Strip, Same Band','Auto-Correlation & Overlay','Next Strip, Same Band';
                   'Prev Strip, Away From Edge','Same Strip, Away from Edge','Next Strip, Away from Edge'};
    subplot(3,3,1);
    plotColors = colormap(hsv(9));
    for j = 1:3
        for k = 1:3
            
            ccActLocalMean(j,k,:) = squeeze(nanmean(nanmean(ccActLocal(:,1:currMaxBands,j,k,:,1),1),2));        
                        
            plotInd = sub2ind([3 3],k,j);
            subplot(3,3,plotInd);      
            hold on
            plot(-p.MaxLag:p.MaxLag,squeeze(ccActLocalMean(j,k,:)),'Color',plotColors(plotInd,:))            
            title(titlesArray{j,k})
            plot(xlim,[0 0],'--k');
            plot([0 0],ylim,'--k')
            
            %Plot them all on the center for comparison
            subplot(3,3,5);
            hold on
            plot(-p.MaxLag:p.MaxLag,squeeze(ccActLocalMean(j,k,:)),'Color',plotColors(plotInd,:))
            
            
        end        
    end
    
    hgsave(ccFigLocal,[p.OutputDirectory filesep 'neighborhood temporal crosscorrelation between acitivities channel ' num2str(p.ChannelIndex(iChan))]);
    
    
    save([p.OutputDirectory filesep 'temporal crosscorrelation channel ' num2str(p.ChannelIndex(iChan))],'ccProtActAll','ccProtActPerBand','ccProtActPerStrip','ccProtActPerWin','ccActLocal','ccActLocalMean')
    
end



