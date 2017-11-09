function [percentageGoodLink,percentageBadLink,percentageMissingLink,cpuTime,tracksFinal,tracksNewFormat] ... 
        = track_stat_3D(movieInfo,tracksSim,tracksSimNF,costMatrices,kalmanFunctions,gapCloseParam,saveFolder,schemeName,frequency)
    
    probDim=3;
    verbose=1;
    saveResults.dir =  [saveFolder filesep schemeName]; %directory where to save input and output
    mkdir(saveResults.dir);
    
    starttime=cputime;    
    [tracksFinal,kalmanInfoLink,errFlag] = ...
        trackCloseGapsKalmanSparse(movieInfo,costMatrices, ... 
                                     gapCloseParam,kalmanFunctions,...
                                     probDim,saveResults,verbose);
    
    

    % computational time
    cpuTime=cputime-starttime;
    
    if(~isempty(tracksFinal))
        tracksNewFormat=TracksHandle(tracksFinal);
        
        % linking percentage
        linkStats=scoreLinksGapsMSLinkType3D(tracksFinal,tracksSim);
        avgperc= nanmean(linkStats(:,3,:)./linkStats(:,1,:))*100;
        percentageGoodLink=avgperc;
        avgperc= nanmean(linkStats(:,4,:)./linkStats(:,1,:))*100;
        percentageBadLink=avgperc;
        avgperc= nanmean( (linkStats(:,1,:)-linkStats(:,3,:))./linkStats(:,1,:))*100;
        percentageMissingLink=avgperc;
        
    else
        tracksNewFormat=[];
        percentageGoodLink=zeros(1,4);
        percentageBadLink=zeros(1,4);
        percentageMissingLink=zeros(1,4);
    end
    
end