function [motionDir,angleWithProt,f2fDisp,paraDirDisp,perpDirDisp,...
    paraProtDisp,perpProtDisp,asymParam] = trackMotionCharProtrusion(...
    tracksFinal,protVecUnit,trackWindowAssign,minLength)

%get number of tracks
numTracks = length(tracksFinal);

%get number of segments making each track and the row of the first
%segment of each track if all track segments were put together
numSegPerTrack = zeros(numTracks,1);
for iTrack = 1 : numTracks
    numSegPerTrack(iTrack) = size(tracksFinal(iTrack).tracksCoordAmpCG,1);
end
trackStartRow = [0; cumsum(numSegPerTrack)];
trackStartRow = 1 + trackStartRow(1:end-1);

%get total number of segments
numSegments = sum(numSegPerTrack);

%reserve memory for output parameters
[f2fDisp,angleWithProt,asymParam] = deal(NaN(numSegments,1));
[motionDir,paraDirDisp,perpDirDisp,paraProtDisp,perpProtDisp] ...
    = deal(NaN(numSegments,2));

%initialize global segment index
iSeg = 0;

%go over all compound tracks
for iTrack = 1 : numTracks
    
    %get current track's coordinates
    trackCoordCurrent = tracksFinal(iTrack).tracksCoordAmpCG;
    xCoord = trackCoordCurrent(:,1:8:end);
    yCoord = trackCoordCurrent(:,2:8:end);
    
    %calculate current track's displacements along x and y
    xCoordDelta = diff(xCoord,[],2);
    yCoordDelta = diff(yCoord,[],2);
    
    %get number of segments in this compound track
    numSeg = numSegPerTrack(iTrack);
    
    %determine which segments are of sufficient length
    segLft = getTrackSEL(trackCoordCurrent);
    indxGood = find(segLft(:,3) >= minLength);
    indxBad  = setdiff(1:numSeg,indxGood);
    
    %calculate average frame-to-frame displacement
    f2fDispCurrent = nanmean( sqrt( xCoordDelta.^2 + yCoordDelta.^2 ) ,2);
    f2fDispCurrent(indxBad) = NaN;
    f2fDisp(iSeg+1:iSeg+numSeg) = f2fDispCurrent;
    
    %reserve memory for this track
    [eigValRatio,angleWithProtCurrent] = deal(NaN(numSeg,1));
    [eigVecMax,paraDirDispCurrent,perpDirDispCurrent,paraProtDispCurrent,...
        perpProtDispCurrent] = deal(NaN(numSeg,2));
    
    %go over segments in track
    for iSegment = indxGood'
        
        %get frame-to-frame displacements along x and y
        trackDisp = [xCoordDelta(iSegment,:); yCoordDelta(iSegment,:)]';
        numEntries = size(trackDisp,1);
        
        %decompose segment's motion to estimate its asymmetry
        posCov = nancov([xCoord(iSegment,:); yCoord(iSegment,:)]');
        [eigVec,eigVal] = eig(posCov);
        eigVal = diag(eigVal);
        eigVal(eigVal<eps) = eps;
        eigValMax = max(eigVal);
        eigValRatio(iSegment) = eigValMax / min(eigVal); %eigenvalue ratio as measure of asymmetry
        
        %get motion direction as eigenvector corresponding to the maximum
        %eigenvalue
        %the direction has to be such that the net displacement along
        %direction of motion is positive
        paraDir = eigVec(:,eigVal==eigValMax)';
        paraTrackDisp = dot(trackDisp,repmat(paraDir,numEntries,1),2); %displacement
        netDispPara = nansum(paraTrackDisp); %net displacement
        eigVecMax(iSegment,:) = sign(netDispPara) * paraDir;
                
        %decompose displacement relative to direction of motion
        paraDir = eigVecMax(iSegment,:);
        perpDir = [-paraDir(2) paraDir(1)];
        paraTrackDisp = dot(trackDisp,repmat(paraDir,numEntries,1),2);
        perpTrackDisp = dot(trackDisp,repmat(perpDir,numEntries,1),2);
        paraDirDispCurrent(iSegment,:) = [nanmean(paraTrackDisp) ...
            nanmean(abs(paraTrackDisp))];
        paraDirDispCurrent(iSegment,paraDirDispCurrent(iSegment,:)==0) = eps;
        %         paraDirDispCurrent(iSegment,:) = max([abs(nanmean(paraTrackDisp)) ...
        %             nanmean(abs(paraTrackDisp))],eps);
        perpDirDispCurrent(iSegment,:) = max([abs(nanmean(perpTrackDisp)) ...
            nanmean(abs(perpTrackDisp))],eps);
                
        %find which window this track falls in
        bigIndx = trackStartRow(iTrack) + iSegment - 1;
        iPara = trackWindowAssign(bigIndx,2);
        iFrame = trackWindowAssign(bigIndx,3);
        
        %if this track falls into some window
        if ~isnan(iPara) && ~isnan(iFrame)
            
            %get protrusion vector
            protParaDir = squeeze(protVecUnit(iPara,iFrame,:))';
            protPerpDir = [-protParaDir(2) protParaDir(1)];
            
            %calculate angle with protrusion vector
            angleWithProtCurrent(iSegment) = acos(dot(eigVecMax(iSegment,:),protParaDir)) * 180 / pi;
            
            %decompose displacement relative to protrusion vector
            paraTrackDisp = dot(trackDisp,repmat(protParaDir,numEntries,1),2);
            perpTrackDisp = dot(trackDisp,repmat(protPerpDir,numEntries,1),2);
            paraProtDispCurrent(iSegment,:) = [nanmean(paraTrackDisp) ...
                nanmean(abs(paraTrackDisp))];
            paraProtDispCurrent(iSegment,paraProtDispCurrent(iSegment,:)==0) = eps;
            %             paraProtDispCurrent(iSegment,:) = max([abs(nanmean(paraTrackDisp)) ...
            %                 nanmean(abs(paraTrackDisp))],eps);
            perpProtDispCurrent(iSegment,:) = max([abs(nanmean(perpTrackDisp)) ...
                nanmean(abs(perpTrackDisp))],eps);
            
        end
        
    end
    
    %store values in output variables
    angleWithProt(iSeg+1:iSeg+numSeg) = angleWithProtCurrent;
    asymParam(iSeg+1:iSeg+numSeg) = eigValRatio;
    motionDir(iSeg+1:iSeg+numSeg,:) = eigVecMax;
    paraDirDisp(iSeg+1:iSeg+numSeg,:) = paraDirDispCurrent;
    perpDirDisp(iSeg+1:iSeg+numSeg,:) = perpDirDispCurrent;
    paraProtDisp(iSeg+1:iSeg+numSeg,:) = paraProtDispCurrent;
    perpProtDisp(iSeg+1:iSeg+numSeg,:) = perpProtDispCurrent;
    
    %update global segment index
    iSeg = iSeg + numSeg;
    
end

disp('')

%% ~~~ the end ~~~
