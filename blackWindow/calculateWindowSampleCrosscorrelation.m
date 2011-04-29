function movieData = calculateWindowSampleCrosscorrelation(movieData,paramIn)


p = paramIn;

%TEMP UNDER CONSTRUCTION
p.ChannelIndex = 3;

%p.BatchMode = false;
p.MaxLag = round(movieData.nFrames_/4);%Gives approximately the maximum possible maxlag
%p.OutputDirectory = [movieData.outputDirectory_ filesep 'sample_autocorrelation_7x3_viscoElastic_MQ10'];
if ~isfield(p,'OutputDirectory')
    error('No output directory specified!')
end

mkClrDir(p.OutputDirectory);

iWSProc = movieData.getProcessIndex('WindowSamplingProcess',1,~p.BatchMode);

if isempty(iWSProc)
    error('The movie must have valid window samples! Please run sampleMovieWindows.m first!');
end

iProtProc = movieData.getProcessIndex('ProtrusionSamplingProcess',1,~p.BatchMode);

if isempty(iProtProc) || ~movieData.processes_{iProtProc}.checkChannelOutput
    disp('Cannot calculate protrusion-cross correlation: No protrusion calculation found!')
else
    
    protSamples = movieData.processes_{iProtProc}.loadChannelOutput;    
    
end
    
nChan = numel(p.ChannelIndex);


for iChan = 1:nChan
    
    %Load the activity samples for this channel
    actSamples = movieData.processes_{iWSProc}.loadChannelOutput(p.ChannelIndex(iChan));
    
    if ndims(actSamples.avg) == 4
        actSamples.avg = squeeze(actSamples.avg);
    else
        error('Finish converting this function to single-object Hunter!')
    end
    
    
    [nStripMax,nBandMax,~] = size(actSamples.avg);

    p.MaxBands = 5;%Maximum number of bands to include in combined cross-correlation analysis. Usually the correlations fall of with distance from edge, so decreasing this number gives stronger correlations.
    
    nObs = zeros(nStripMax,nBandMax);
    %TEMP - this initialization doesn't seem to convince mLint - whats
    %wrong????? HLE
    actTimeSeries(1:nStripMax,1:nBandMax) = struct('observations',[]);
    protTimeSeries(1:nStripMax) = struct('observations',[]);
    ccProtActPerWin = nan(nStripMax,nBandMax,2*p.MaxLag+1,2);
    ccProtActPerStrip = nan(nStripMax,2*p.MaxLag+1,2);
    ccProtActPerBand = nan(nBandMax,2*p.MaxLag+1,2);
    
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
            
            if (nObs(k,l) - p.MaxLag) >= 3*p.MaxLag                                 
                ccProtActPerWin(k,l,:,:) = crossCorr(actTimeSeries(k,l),protTimeSeries(k,l),p.MaxLag);                                                          
                
            end                                   
            
        end
        
        if (sum(squeeze(nObs(k,nObs(k,1:p.MaxBands)>p.MaxLag))) - p.MaxLag) >= (3*p.MaxLag) %Make sure we have enough points
            ccProtActPerStrip(k,:,:) = crossCorr(actTimeSeries(k,1:p.MaxBands),protTimeSeries(k,1:p.MaxBands),p.MaxLag);%we need at least 3 points per traj for trend removal
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
    
    
    ccProtActAll = crossCorr(actTimeSeries(:,1:p.MaxBands),protTimeSeries,p.MaxLag);       
    
    if p.BatchMode
        ccFig = figure('Visible','off');
    else
        ccFig = figure;
    end

    plot(-p.MaxLag:p.MaxLag,ccProtActAll(:,1),'r');
    hold on    
    plotTransparent(-p.MaxLag:p.MaxLag,ccProtActAll(:,1),ccProtActAll(:,2),'r',.2,0);
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
            
            ccActLocalMean(j,k,:) = squeeze(nanmean(nanmean(ccActLocal(:,1:p.MaxBands,j,k,:,1),1),2));        
                        
            plotInd = sub2ind([3 3],k,j);
            subplot(3,3,plotInd);      
            hold on
            plot(-p.MaxLag:p.MaxLag,squeeze(ccActLocalMean(j,k,:)),'Color',plotColors(plotInd,:))            
            title(titlesArray{j,k})
            plot(xlim,[0 0],'--k');
            
            %Plot them all on the center for comparison
            subplot(3,3,5);
            hold on
            plot(-p.MaxLag:p.MaxLag,squeeze(ccActLocalMean(j,k,:)),'Color',plotColors(plotInd,:))
            
            
        end        
    end
    
    hgsave(ccFigLocal,[p.OutputDirectory filesep 'neighborhood temporal crosscorrelation between acitivities channel ' num2str(p.ChannelIndex(iChan))]);
    
    
    save([p.OutputDirectory filesep 'temporal crosscorrelation channel ' num2str(p.ChannelIndex(iChan))],'ccProtActAll','ccProtActPerBand','ccProtActPerStrip','ccProtActPerWin','ccActLocal','ccActLocalMean')
    
end



