function movieData = calculateWindowSampleAutocorrelation(movieData,paramIn)


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

nChan = numel(p.ChannelIndex);



for iChan = 1:nChan

    samples = movieData.processes_{iWSProc}.loadChannelOutput(p.ChannelIndex(iChan));
    
    [nObjMax,nStripMax,nBandMax,nFrames] = size(samples.avg);
    
    actTimeSeries(1:nObjMax,1:nStripMax,1:nBandMax) = struct('observations',[]);
    nObs = nan(nObjMax,nStripMax,nBandMax);
    acPerWin = nan(nObjMax,nStripMax,nBandMax,p.MaxLag+1,2);
    acPerBand = nan(nObjMax,p.MaxLag+1,2,nBandMax);
    acPerStrip = nan(nObjMax,p.MaxLag+1,2,nStripMax);
    
    for j = 1:nObjMax
        
        for l = 1:nBandMax
            
            for k = 1:nStripMax
                %Extract the current time-series
                actTimeSeries(j,k,l).observations = squeeze(samples.avg(j,k,l,:));
                
                %Leave the empty series as placeholders
                nObs(j,k,l) = nnz(~isnan(actTimeSeries(j,k,l).observations));                
                
                %De-trend the current time-series.
                if nObs(j,k,l) >= 3 %We need at least 3 points for de-trending to make sense                    
                    actTimeSeries(j,k,l).observations = removeLinearTrend(actTimeSeries(j,k,l).observations,[],0);
                    acPerWin(j,k,l,:,:) = autoCorr(actTimeSeries(j,k,l),p.MaxLag,-1);                    
                end
                
            end
            if (sum(squeeze(nObs(j,nObs(j,:,l)>p.MaxLag,l))) - p.MaxLag) >= (3*p.MaxLag) %Make sure we have enough points
                acPerBand(j,:,:,l) = autoCorr(actTimeSeries(j,nObs(j,:,l)>=3,l),p.MaxLag,-1);%we need at least 3 points per traj for trend removal
            end
            
        end                
        %I think this is a dumb way to do this but I'm too tired to tell.
        for k = 1:nStripMax
            if (sum(squeeze(nObs(j,k,nObs(j,k,:)>p.MaxLag))) - p.MaxLag) >= (3*p.MaxLag) %Make sure we have enough points by K's criteria
                acPerStrip(j,:,:,k) = autoCorr(actTimeSeries(j,k,nObs(j,k,:)>=3),p.MaxLag,-1);
            end
        end
        
    end
    
    nPtsTot = sum(nObs(:));
    
    [acActTimeAll,err(1)] = autoCorr(actTimeSeries(nObs(:)>3),p.MaxLag,-1);

    if p.BatchMode
        acFig = figure('Visible','off');
    else
        acFig = figure;
    end

    plot(0:p.MaxLag,acActTimeAll(:,1),'r');
    hold on
    plot(xlim,ones(1,2)*1.96/sqrt(nPtsTot),'--k');
    plot(xlim,-ones(1,2)*1.96/sqrt(nPtsTot),'--k');
    plotTransparent(0:p.MaxLag,acActTimeAll(:,1),acActTimeAll(:,2),'r',.2,0);
    xlabel('Time Lag, frames'); %TEMP - get time interval and scale if available!!!!
    ylabel('Autocorrelation');
    title(['All Samples Combined, Temporal Autocorrelation, Channel ' num2str(p.ChannelIndex(iChan))]);
    if any(err)
        error('Problem calculating autocorrelation!')
    end

    hgsave(acFig,[p.OutputDirectory filesep 'combined temporal autocorrelation channel ' num2str(p.ChannelIndex(iChan))]);            
    
    if p.BatchMode
        acFigBand = figure('Visible','off');
    else
        acFigBand = figure;
    end
    
    
    
    surf(squeeze(acPerBand(1,:,1,:)))%TEMP - needs to be general to n Objects>1
    view([-125 20])
    xlabel('Into Cell (band #)')
    ylabel('Time Lag (frames)')
    zlabel('Autocorrelation')
    
    title(['Per Band, Temporal Autocorrelation, Channel ' num2str(p.ChannelIndex(iChan))]);
    
    hgsave(acFigBand,[p.OutputDirectory filesep 'per band temporal autocorrelation channel ' num2str(p.ChannelIndex(iChan))]);
    
    if p.BatchMode
        acFigStrip = figure('Visible','off');
    else
        acFigStrip = figure;
    end
            
    
    surf(squeeze(acPerStrip(1,:,1,:)))%TEMP - needs to be general to n Objects>1
    view([-125 20])
    xlabel('Along Cell Edge (strip #)')
    ylabel('Time Lag (frames)')
    zlabel('Autocorrelation')
    
    title(['Per Strip, Temporal Autocorrelation, Channel ' num2str(p.ChannelIndex(iChan))]);
    
    hgsave(acFigStrip,[p.OutputDirectory filesep 'per strip temporal autocorrelation channel ' num2str(p.ChannelIndex(iChan))]);
        
    
    if p.BatchMode
        acFigWin = figure('Visible','off');
    else
        acFigWin = figure;
    end
    acSum = sqrt(nansum(acPerWin(:,:,:,:,1) .^2,4));
    nCol = 1024;
    winColors = colormap(jet(nCol));
    maxAc = mean(acSum(acSum(:)>0))+std(acSum(acSum(:)>0));
    minAc = 1;
    for j = 1:nObjMax
        for k = 1:nStripMax
            for l = 1:nBandMax
                
                if nObs(j,k,l) >= 3
                    iCol = round(interp1(linspace(minAc,maxAc,nCol),1:nCol,min(acSum(j,k,l),maxAc)));%This is probably slow - do it manually?                    
                    plotString{j}{k}{l} = {winColors(iCol,:)};
                else
                    plotString{j}{k}{l} = {'k','FaceAlpha',0};
                end
                
            end
        end
    end
    
    iWinProc = movieData.getProcessIndex('WindowingProcess',1,0);
    windows = movieData.processes_{iWinProc}.loadChannelOutput(1);
    
    plotWindows(windows,plotString)
    axis off
    caxis([minAc maxAc])
    colorbar
    title('Total Autocorrelation at all Lags, Per window')
    
    hgsave(acFigWin,[p.OutputDirectory filesep 'per window temporal autocorrelation channel ' num2str(p.ChannelIndex(iChan))]);
            
    save([p.OutputDirectory filesep 'temporal autocorrelation channel ' num2str(p.ChannelIndex(iChan))],'acActTimeAll','acPerBand','acPerStrip','acPerWin','nObs');   
    
end

    



