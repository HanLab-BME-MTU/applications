function [projDataBatch]=metaEB3analysis(projList,secPerFrame,pixSizeNm)


for i=1:length(projList)
    
    if exist([projList(i,1).anDir filesep 'testTrack.mat'])~=0
        
        projData.imDir=projList(i,1).imDir;
        projData.anDir=projList(i,1).anDir;
        projData.secPerFrame=secPerFrame;
        projData.pixSizeNm=pixSizeNm;
       
        
        load([projList(i,1).anDir filesep 'testTrack.mat'])
        nTracks=length(tracksFinal);
        projData.nTracks=nTracks;
        [trackedFeatureInfoInterp,trackInfo,trackVelocities] = getVelocitiesFromMat(tracksFinal,3);
        
        [projData.trackVelocities] = pixPerFrame2umPerMin(trackVelocities.segmentAvgs,secPerFrame,pixSizeNm);
        projData.xCoord=trackedFeatureInfoInterp(:,1:8:end);
        projData.yCoord=trackedFeatureInfoInterp(:,2:8:end);
        
        %concatenate all segments and gaps into n x 4 matrices, then add...
        %each has format [trackNum startFrame endFrame velocity seg/gapType trackLengthFrames]
        segs   = vertcat(trackInfo.seg);  segs =  [segs   1*ones(size(segs,1),1)   segs(:,3)-segs(:,2)];
        fgaps  = vertcat(trackInfo.fgap); fgaps = [fgaps  2*ones(size(fgaps,1),1) fgaps(:,3)-fgaps(:,2)];
        bgaps  = vertcat(trackInfo.bgap); bgaps = [bgaps  3*ones(size(bgaps,1),1) bgaps(:,3)-bgaps(:,2)];
        ugaps  = vertcat(trackInfo.ugap); ugaps = [ugaps  4*ones(size(ugaps,1),1) ugaps(:,3)-ugaps(:,2)];
        
        projData.meanSegVel = mean(pixPerFrame2umPerMin(segs(:,4),secPerFrame,pixSizeNm));
        projData.medianSegVel = median(pixPerFrame2umPerMin(segs(:,4),secPerFrame,pixSizeNm));        
        projData.perc95SegVel = prctile(pixPerFrame2umPerMin(segs(:,4),secPerFrame,pixSizeNm),95);

        % sort to see track profiles in order
        aT=sortrows([segs; fgaps; bgaps; ugaps],[1 2]);
        
        % convert pix/frame --> micron/min velocities
        [aT(:,4)] = pixPerFrame2umPerMin(aT(:,4),secPerFrame,pixSizeNm);
        
        tS=zeros(nTracks,7); %temp for trackStatistics
        for j=1:nTracks
            tS(j,1)=j;
            idx=find(aT(:,1)==j);
            tS(j,2)=max(aT(idx,3))-min(aT(idx,2)); % total track length in frames
            tS(j,3)=length(idx); % number of segs and gaps that make up the track
            
            tS(j,4)=sum(aT(idx(aT(idx,5)==1),6)); % total seg frames
            tS(j,5)=sum(aT(idx(aT(idx,5)==2),6)); % total seg frames
            tS(j,6)=sum(aT(idx(aT(idx,5)==3),6)); % total seg frames
            tS(j,7)=sum(aT(idx(aT(idx,5)==4),6)); % total seg frames
        
        end

        projData.allTracks=aT;
        projData.trackStats=tS;
        
        % save each projData in its own directory
        save([projList(i,1).anDir filesep 'projData'],'projData')

        % keep a structure for all the projects
        projDataBatch(i,1)=projData;
        
    else
        disp([projList(i,1).anDir ' has not been analyzed'])
    end
end


popVelScatter(projDataBatch)
