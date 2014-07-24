function postProcessMovieArrayObjectTracking(MA)


%UNDER CONSTRUCTION YOU ASSHOLES!!!

nMovies = numel(MA);

p.BatchMode = false;
p.ChannelIndex = 1;
p.SmoothPar = 1e-9;


for iMov = 1:nMovies
    
    
    iProc = MA(iMov).getProcessIndex('MaskObjectTrackingProcess',1,0);    
    
    
    if ~isempty(iProc)
        
        tData{iMov} = 0:MA(iMov).timeInterval_:(MA(iMov).timeInterval_ * (MA(iMov).nFrames_-1));
    
        %Already scaled in Z, so we just multiply by XY pixel size.
        objTracks{iMov} = squeeze(MA(iMov).processes_{iProc}.loadChannelOutput(p.ChannelIndex)) * MA(iMov).pixelSize_;
        objTrackSpline(iMov) = csaps(tData{iMov},objTracks{iMov}',p.SmoothPar);
        objTrackVel(iMov) = fnder(objTrackSpline(iMov),1);
        
        if p.BatchMode
            pmTrackFig = figure('Visible','off');
        else
            pmTrackFig = figure;
        end
        hold on;
        plot3(objTracks{iMov}(:,1),objTracks{iMov}(:,2),objTracks{iMov}(:,3),'r')
        axis equal
        fnplt(objTrackSpline(iMov))
        
        if p.BatchMode
            pmVelFig = figure('Visible','off');
        else
            pmVelFig = figure;
        end
        subplot(2,1,1)
        hold on;
        spVals = ppval(objTrackSpline(iMov),tData{iMov});        
        ylabel('Displacement (nm)')
        xlabel('Time (s)')
        plot(tData{iMov},objTracks{iMov}(:,1) - objTracks{iMov}(1,1))                
        plot(tData{iMov},objTracks{iMov}(:,2) - objTracks{iMov}(1,2),'r')                
        plot(tData{iMov},objTracks{iMov}(:,3) - objTracks{iMov}(1,3),'g')        
        legend('x','y','z')
        plot(tData{iMov},spVals(1,:)-spVals(1,1),'--b')
        plot(tData{iMov},spVals(2,:)-spVals(2,1),'--r')
        plot(tData{iMov},spVals(3,:)-spVals(3,1),'--g')        
        
        subplot(2,1,2)
        hold on;
        spVals = ppval(objTrackVel(iMov),tData{iMov});        
        ylabel('Velocity (nm/s)')
        xlabel('Time (s)')
        plot(tData{iMov}(1:end-1),diff(objTracks{iMov}(:,1)))                
        plot(tData{iMov}(1:end-1),diff(objTracks{iMov}(:,2)),'r')
        plot(tData{iMov}(1:end-1),diff(objTracks{iMov}(:,3)),'g')
        legend('x','y','z')
        plot(tData{iMov},spVals(1,:)-spVals(1,1),'--b')
        plot(tData{iMov},spVals(2,:)-spVals(2,1),'--r')
        plot(tData{iMov},spVals(3,:)-spVals(3,1),'--g')                                
        
        
    else        
        disp(['Skipping movie ' num2str(iMov) ' : No object tracking found!'])
    end
    
    
    
    
end



