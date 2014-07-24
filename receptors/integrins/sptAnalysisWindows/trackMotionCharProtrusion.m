function trackChar = trackMotionCharProtrusion(tracksFinal,protSamples,...
    trackWindowAssignComp,minLength)
%TRACKMOTIONCHARPROTRUSION calculate various trajectory motion characteristics, also relative to cell protrusion vector
%
%SYNOPSIS trackChar = trackMotionCharProtrusion(tracksFinal,protSamples,...
%    trackWindowAssignComp,minLength)
%
%INPUT  tracksFinal    : The tracks as output by trackCloseGapsKalman.
%       protSamples    : The protrusion samples as output by the windowing
%                        software.
%                        Optional. If not input, characteristics relative
%                        to protrusion vector are not calculated.
%       trackWindowAssignComp: Structure array indicating the window
%                        assignment of each track.
%                        Optional. If not input, characteristics relative
%                        to protrusion vector are not calculated.
%       minLength      : Minimum length of a trajectory to include in
%                        analysis.
%                        Optional. Default: 5.
%
%OUTPUT trackChar      : Structure array of the same length as tracksFinal.
%                        Contains the following fields:
%           .motionDir     : Track's direction of motion.
%           .angleWithProt : Angle between track's direction of motion and
%                            cell protrusion vector.
%           .f2fDisp       : Average frame-to-frame displacement.
%           .paraDirDisp   : Average displacement parallel to track
%                            direction of motion.
%           .perpDirDisp   : Average displacement perpendicular to track
%                            direction of motion.
%           .paraProtDisp  : Average displacement parallel to cell
%                            protrusion vector.
%           .perpProtDisp  : Average displacement perpendicult to cell
%                            protrusion vector.
%           .asymParam     : Measure of track asymmetry.
%           .f2fDispRMS    : Root-mean-square frame-to-frame displacement.
%
%Khuloud Jaqaman, January 2012

%% Input

if nargin < 2 || isempty(protSamples) || nargin < 3 || isempty(trackWindowAssignComp)
    charProt = 0;
else
    charProt = 1;
end

if nargin < 4 || isempty(minLength)
    minLength = 5;
end

%calculate unit protrusion vectors
protVec = protSamples.avgVector;
protVecMag = sqrt(sum(protVec.^2,3));
protVecUnit = protVec ./ repmat(protVecMag,[1 1 2]);

%% Track characteristics

%get number of tracks
numTracks = length(tracksFinal);

%reserve memory for output parameters
trackChar = repmat(struct('motionDir',[],'angleWithProt',[],'f2fDisp',[],...
    'paraDirDisp',[],'perpDirDisp',[],'paraProtDisp',[],'perpProtDisp',[],...
    'asymParam',[],'f2fDispRMS',[]),numTracks,1);

%go over all compound tracks
parfor iTrack = 1 : numTracks

    %get current track's coordinates
    trackCoordCurrent = tracksFinal(iTrack).tracksCoordAmpCG;
    xCoord = trackCoordCurrent(:,1:8:end);
    yCoord = trackCoordCurrent(:,2:8:end);

    %calculate current track's displacements along x and y
    xCoordDelta = diff(xCoord,[],2);
    yCoordDelta = diff(yCoord,[],2);

    %get number of segments in this compound track
    numSeg = size(xCoord,1);

    %determine which segments are of sufficient length
    segLft = getTrackSEL(trackCoordCurrent);
    indxGood = find(segLft(:,3) >= minLength);
    indxBad  = setdiff(1:numSeg,indxGood);

    %calculate average frame-to-frame displacement
    f2fDispCurrent = nanmean( sqrt( xCoordDelta.^2 + yCoordDelta.^2 ) ,2);
    f2fDispCurrent(indxBad) = NaN;
    trackChar(iTrack).f2fDisp = f2fDispCurrent;

    %calculate root mean square frame-to-frame displacement
    f2fDispRMSCurrent = sqrt( nanmean(xCoordDelta.^2+yCoordDelta.^2,2) );
    f2fDispRMSCurrent(indxBad) = NaN;
    trackChar(iTrack).f2fDispRMS = f2fDispRMSCurrent;

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

        if charProt
            
            %find which window this track falls in
            iPara = trackWindowAssignComp(iTrack).assignment(iSegment,2);
            iFrame = trackWindowAssignComp(iTrack).assignment(iSegment,3);
            
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
                perpProtDispCurrent(iSegment,:) = [abs(nanmean(perpTrackDisp)) ...
                    nanmean(abs(perpTrackDisp))];
                perpProtDispCurrent(iSegment,perpProtDispCurrent(iSegment,:)==0) = eps;
                
            end
            
        end

    end

    %store values in output variables
    trackChar(iTrack).angleWithProt = angleWithProtCurrent;
    trackChar(iTrack).asymParam = eigValRatio;
    trackChar(iTrack).motionDir = eigVecMax;
    trackChar(iTrack).paraDirDisp = paraDirDispCurrent;
    trackChar(iTrack).perpDirDisp = perpDirDispCurrent;
    trackChar(iTrack).paraProtDisp = paraProtDispCurrent;
    trackChar(iTrack).perpProtDisp = perpProtDispCurrent;

end

%% ~~~ the end ~~~
