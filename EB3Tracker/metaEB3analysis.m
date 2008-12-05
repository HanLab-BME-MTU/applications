function [projData]=metaEB3analysis(projList,secPerFrame,pixSizeNm)

projData(length(projList),1).imDir=[];
projData(length(projList),1).anDir=[];

counter=0;
for i=1:length(projList)
    
    if exist([projList(i,1).anDir filesep 'testTrack.mat'])~=0
        counter=counter+1;
        
        projData(i,1).imDir=projList(i,1).imDir;
        projData(i,1).anDir=projList(i,1).anDir;
        projData(i,1).secPerFrame=secPerFrame;
        projData(i,1).pixSizeNm=pixSizeNm;
       
        
        load([projList(i,1).anDir filesep 'testTrack.mat'])
        projData(counter).nTracks=length(tracksFinal);
        [trackedFeatureInfoInterp,trackInfo,trackVelocities] = getVelocitiesFromMat(tracksFinal,3);
        
        [projData(i,1).trackVelocities] = pixPerFrame2umPerMin(trackVelocities.segmentAvgs,secPerFrame,pixSizeNm);
        projData(i,1).xCoord=trackedFeatureInfoInterp(:,1:8:end);
        projData(i,1).yCoord=trackedFeatureInfoInterp(:,2:8:end);
        
        %concatenate all segments and gaps into n x 4 matrices, then add...
        %each has format [trackNum startFrame endFrame velocity seg/gapType trackLengthFrames]
        segs   = vertcat(trackInfo.seg);  segs =  [segs   1*ones(size(segs,1),1)   segs(:,3)-segs(:,2)];
        fgaps  = vertcat(trackInfo.fgap); fgaps = [fgaps  2*ones(size(fgaps,1),1) fgaps(:,3)-fgaps(:,2)];
        bgaps  = vertcat(trackInfo.bgap); bgaps = [bgaps  3*ones(size(bgaps,1),1) bgaps(:,3)-bgaps(:,2)];
        ugaps  = vertcat(trackInfo.ugap); ugaps = [ugaps  4*ones(size(ugaps,1),1) ugaps(:,3)-ugaps(:,2)];
        
        % sort to see track profiles in order
        projData(i,1).allTracks=sortrows([segs; fgaps; bgaps; ugaps],[1 2]);
        
        % convert pix/frame --> micron/min velocities
        [projData(i,1).allTracks(:,4)] = pixPerFrame2umPerMin(projData(i,1).allTracks(:,4),secPerFrame,pixSizeNm);
        temp=projData;
        projData=projData(i,1);
        save([projList(i,1).anDir filesep 'projData'],'projData')
        projData=temp;
        clear temp;
        
    else
        disp([projList(i,1).anDir ' has not been analyzed'])
    end
end

if counter~=length(projList)
    projData=projData(1:counter);
end

popVelScatter(projData)
