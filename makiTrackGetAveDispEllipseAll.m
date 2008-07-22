function [longVecS,longVecE,shortVecS,shortVecE,shortVecS3D,shortVecE3D] = ...
    makiTrackGetAveDispEllipseAll(xyzVel,brownStd,trackType,undetBrownStd,timeWindow,...
    brownStdMult,linStdMult,timeReachConfB,timeReachConfL,minSearchRadius,...
    maxSearchRadius,useLocalDensity,closestDistScale,maxStdMult,...
    nnDistLinkedFeat,nnWindow,trackStartTime,trackEndTime,probDim,...
    trackStartKinType,trackEndKinType)
%MAKITRACKGETAVEDISPELLIPSEALL determines the search ellipse and expected displacement for HeLa kinetochores
%
%SYNOPSIS [longVecS,longVecE,shortVecS,shortVecE,shortVecS3D,shortVecE3D] = ...
%    makiTrackGetAveDispEllipseAll(xyzVel,brownStd,trackType,undetBrownStd,timeWindow,...
%    brownStdMult,linStdMult,timeReachConfB,timeReachConfL,minSearchRadius,...
%    maxSearchRadius,useLocalDensity,closestDistScale,maxStdMult,...
%    nnDistLinkedFeat,nnWindow,trackStartTime,trackEndTime,probDim,...
%    trackStartKinType,trackEndKinType)
%
%INPUT  xyzVel         : Velocity in x, y and z (if 3D).
%       brownStd       : Standard deviation of Brownian motion steps.
%       trackType      : Type of track. 1 for directed, 0 for Brownian, NaN for undetermined.
%       undetBrownStd  : Standard deviation of Brownian motion steps to be used
%                        for undetermined tracks.
%       timeWindow     : Maximum gap size.
%       brownStdMult   : Multiplication factor to go from average Brownian
%                        displacement to search radius.
%       linStdMult     : Multiplication factor to go from average linear
%                        displacement to search radius.
%       timeReachConfB : Time gap for Brownian motion to reach confinement.
%       timeReachConfL : Time gap for linear motion to reach confinement.
%       minSearchRadius: Minimum allowed search radius.
%       maxSearchRadius: Maximum allowed search radius for linking between
%                        two consecutive frames. It will be expanded for
%                        different gap lengths based on the time scaling of
%                        Brownian motion.
%       useLocalDensity: 1 if local density of features is used to expand 
%                        their search radius if possible, 0 otherwise.
%       closestDistScale:Scaling factor of nearest neighbor distance.
%       maxStdMult     : Maximum value of factor multiplying std to get
%                        search radius.
%       nnDistLinkedFeat:Matrix indicating the nearest neighbor
%                        distances of features linked together within
%                        tracks.
%       nnWindow       : Time window to be used in estimating the
%                        nearest neighbor distance of a track at its start
%                        and end.
%       trackStartTime : Starting time of all tracks.
%       trackEndTime   : Ending time of all tracks.
%       probDim        : Problem dimensionality. 2 (for 2D) or 3 (for 3D).
%       trackStartKinType: Type of kinetochore at track start.
%       trackEndKinType: Type of kinetochore at track end.
%
%OUTPUT longVecS  : Vector defining long radius of search ellipse/ellipsoid at the
%                   starts of tracks.
%       longVecE  : Vector defining long radius of search ellipse/ellipsoid at the
%                   ends of tracks.
%       shortVecS : Vector defining short radius of search ellipse/ellipsoid at the
%                   starts of tracks.
%       shortvecE : Vector defining short radius of search ellipse/ellipsoid at the
%                   ends of tracks.
%       shortVecS3D:Vector defining 2nd short radius of search ellipse/ellipsoid at the
%                   starts of tracks in case of 3D.
%       shortvecE3D:Vector defining 2nd short radius of search ellipse/ellipsoid at the
%                   ends of tracks in case of 3D.
%       errFlag   : 0 if function executes normally, 1 otherwise
%
%       dispDrift : Vector of expected displacement along x and y due to
%                   drift (not output any more - although still
%                   calculated).
%       dispBrown : Expected displacement along x (= along y) due to
%                   Brownian motion (not output any more - although still
%                   calculated).
%
%REMARKS Drift is assumed to look more like 1D diffusion, i.e. the particle
%goes back and forth along a line
%
%Khuloud Jaqaman, April 2007

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% dispDrift = [];
% dispBrown = [];
longVecS  = [];
longVecE  = [];
shortVecS = [];
shortVecE = [];
shortVecS3D = [];
shortVecE3D = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check whether correct number of input arguments was used
if nargin ~= nargin('makiTrackGetAveDispEllipseAll')
    disp('--makiTrackGetAveDispEllipseAll: Incorrect number of input arguments!');
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Determine expected displacement and search ellipse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%determine number of tracks
numTracks = size(xyzVel,1);

%reserve memory for output
% dispDrift = zeros(probDim,timeWindow,numTracks);
% dispBrown = zeros(timeWindow,numTracks);
longVecS  = zeros(probDim,timeWindow,numTracks);
longVecE  = zeros(probDim,timeWindow,numTracks);
shortVecS = zeros(probDim,timeWindow,numTracks);
shortVecE = zeros(probDim,timeWindow,numTracks);

%define square root of "problem dimension" to avoid calculating it many times
sqrtDim = sqrt(probDim);

%put time scaling of linear motion in a vector
timeScalingLin = [sqrt(1:timeReachConfL) sqrt(timeReachConfL) * ...
    (2:timeWindow-timeReachConfL+1).^0.01];

%put time scaling of Brownian motion in a vector
timeScalingBrown = [sqrt(1:timeReachConfB) sqrt(timeReachConfB) * ...
    (2:timeWindow-timeReachConfB+1).^0.01];

%scale maxSearchRadius like Brownian motion (it's only imposed on the
%Brownian aspect of tracks)
maxSearchRadius = repmat(maxSearchRadius',1,timeWindow) .* repmat(timeScalingBrown,3,1);

%determine the nearest neighbor distances of tracks at their starts and ends
windowLimS = min([trackStartTime+nnWindow trackEndTime],[],2);
windowLimE = max([trackEndTime-nnWindow trackStartTime],[],2);
nnDistTracksS = zeros(numTracks,1);
nnDistTracksE = zeros(numTracks,1);
for iTrack = 1 : numTracks
    nnDistTracksS(iTrack) = min(nnDistLinkedFeat(iTrack,...
        trackStartTime(iTrack):windowLimS(iTrack)));
    nnDistTracksE(iTrack) = min(nnDistLinkedFeat(iTrack,...
        windowLimE(iTrack):trackEndTime(iTrack)));
end

for iTrack = 1 : numTracks
    
    %get kinetochore types at start and end of track
    kinTypeStart = trackStartKinType(iTrack);
    kinTypeEnd = trackEndKinType(iTrack);
    
    switch trackType(iTrack)

        case 1
            
            %get velocity, its magnitude and direction of motion
            velDrift = xyzVel(iTrack,:)';
            velMag = sqrt(velDrift' * velDrift);
            directionMotion = velDrift / velMag;
            
            %obtain vector(s) perpendicular to direction of motion
            if probDim == 2 %in 2D case, 1 vector needed
                perpendicular = [-directionMotion(2) directionMotion(1)]';
            else %in 3D case, 2 vectors needed
                perpendicular = [-directionMotion(2) directionMotion(1) 0]';
                perpendicular = perpendicular / (sqrt(perpendicular'*perpendicular));
                perpendicular3D = cross(directionMotion,perpendicular);
            end

            %calculate the expected displacement due to drift for all time
            %gaps
            dispDrift1 = repmat(velDrift,1,timeWindow) .* repmat(timeScalingLin,probDim,1);

            %calculate the expected displacement along x (= along y, [z]) due to
            %brownian motion for all time gaps
            dispBrown1 = brownStd(iTrack) * timeScalingBrown;
            
            %copy brownStdMult into vector that might be modified using
            %local density
            brownStdMultModS = brownStdMult'; %for track start
            brownStdMultModE = brownStdMult'; %for track end

            %if local density information is used to expand search radius ...
            if useLocalDensity

                %divide the track's nearest neighbor distance at its start
                %/closestDistScale by expected Brownian displacement
                ratioDist2Std = repmat(nnDistTracksS(iTrack)/closestDistScale,...
                    1,timeWindow) ./ dispBrown1;

                %make ratios larger than maxStdMult equal to maxStdMult
                ratioDist2Std(ratioDist2Std > maxStdMult) = maxStdMult;

                %expand search radius multiplication factor at track start
                %if possible
                brownStdMultModS = max([brownStdMultModS; ratioDist2Std]);

                %divide the track's nearest neighbor distance at its end
                %/closestDistScale by expected Brownian displacement
                ratioDist2Std = repmat(nnDistTracksE(iTrack)/closestDistScale,...
                    1,timeWindow) ./ dispBrown1;

                %make ratios larger than maxStdMult equal to maxStdMult
                ratioDist2Std(ratioDist2Std > maxStdMult) = maxStdMult;

                %expand search radius multiplication factor at track end
                %if possible
                brownStdMultModE = max([brownStdMultModE; ratioDist2Std]);

            end

            %determine the long vectors of the search ellipses for all time
            %gaps
            longVec1 = repmat(linStdMult',probDim,1) .* dispDrift1 + ...
                repmat((brownStdMult' .* dispBrown1 * sqrtDim),probDim,1) .* ...
                repmat(directionMotion,1,timeWindow);
            longVecMag = sqrt((diag(longVec1' * longVec1))');  %magnitude
            longVecDir = longVec1 ./ repmat(longVecMag,probDim,1); %direction
            
            %determine the short vectors at track starts
            shortVecS1 = repmat((brownStdMultModS .* dispBrown1 * sqrtDim),probDim,1) .* ...
                repmat(perpendicular,1,timeWindow);
            shortVecSMag = sqrt((diag(shortVecS1' * shortVecS1))');  %magnitude
            shortVecSDir = shortVecS1 ./ repmat(shortVecSMag,probDim,1); %direction

            %determine the short vectors at track ends
            shortVecE1 = repmat((brownStdMultModE .* dispBrown1 * sqrtDim),probDim,1) .* ...
                repmat(perpendicular,1,timeWindow);
            shortVecEMag = sqrt((diag(shortVecE1' * shortVecE1))');  %magnitude
            shortVecEDir = shortVecE1 ./ repmat(shortVecEMag,probDim,1); %direction
            
            %             %output the absolute value of dispDrift1
            %             dispDrift1 = abs(dispDrift1);

            %make sure that long vectors are longer than minimum allowed
            %starts:
            longVecSMag = max([longVecMag;repmat(minSearchRadius(kinTypeStart+1),1,timeWindow)]); %compare to minimum
            longVecS1 = repmat(longVecSMag,probDim,1) .* longVecDir; %new long vector
            longVecEMag = max([longVecMag;repmat(minSearchRadius(kinTypeEnd+1),1,timeWindow)]); %compare to minimum
            longVecE1 = repmat(longVecEMag,probDim,1) .* longVecDir; %new long vector

            %make sure that short vectors at track starts are within
            %allowed range
            shortVecSMag = max([shortVecSMag;repmat(minSearchRadius(kinTypeStart+1),1,timeWindow)]); %compare to minimum
            shortVecSMag = min([shortVecSMag;maxSearchRadius(kinTypeStart+1,:)]); %compare to maximum
            shortVecS1 = repmat(shortVecSMag,probDim,1) .* shortVecSDir; %new short vector

            %make sure that short vectors at track ends are within allowed
            %range
            shortVecEMag = max([shortVecEMag;repmat(minSearchRadius(kinTypeEnd+1),1,timeWindow)]); %compare to minimum
            shortVecEMag = min([shortVecEMag;maxSearchRadius(kinTypeEnd+1,:)]); %compare to maximum
            shortVecE1 = repmat(shortVecEMag,probDim,1) .* shortVecEDir; %new short vector
            
            %save values for this track
            %             dispDrift(:,:,iTrack) = dispDrift1;
            %             dispBrown(:,iTrack) = dispBrown1;
            longVecS(:,:,iTrack) = longVecS1;
            longVecE(:,:,iTrack) = longVecE1;
            shortVecS(:,:,iTrack) = shortVecS1;
            shortVecE(:,:,iTrack) = shortVecE1;
            
            %construct additional short vectors for 3D problems
            if probDim == 3
                shortVecS13D = repmat(shortVecSMag,probDim,1) .* repmat(perpendicular3D,1,timeWindow);
                shortVecE13D = repmat(shortVecEMag,probDim,1) .* repmat(perpendicular3D,1,timeWindow);
                shortVecS3D(:,:,iTrack) = shortVecS13D;
                shortVecE3D(:,:,iTrack) = shortVecE13D;
            end

        case 0
            
            %take direction of motion to be along x and construct
            %perpendicular(s)
            if probDim == 2
                directionMotion = [1 0]';
                perpendicular = [0 1]';
            else
                directionMotion = [1 0 0]';
                perpendicular = [0 1 0]';
                perpendicular3D = [0 0 1]';
            end

            %calculate the expected displacement along x (= along y, [z]) due to
            %brownian motion for all time gaps
            dispBrown1 = brownStd(iTrack) * timeScalingBrown;

            %copy brownStdMult into vector that might be modified using
            %local density
            brownStdMultModS = brownStdMult'; %for track start
            brownStdMultModE = brownStdMult'; %for track end
            
            %if local density information is used to expand search radius ...
            if useLocalDensity

                %divide the track's nearest neighbor distance/closestDistScale
                %by expected Brownian displacement
                ratioDist2Std = repmat(nnDistTracksS(iTrack)/closestDistScale,...
                    1,timeWindow) ./ dispBrown1;

                %make ratios larger than maxStdMult equal to maxStdMult
                ratioDist2Std(ratioDist2Std > maxStdMult) = maxStdMult;

                %expand start's search radius multiplication factor if possible
                brownStdMultModS = max([brownStdMultModS; ratioDist2Std]);

                %divide the track's nearest neighbor distance at its end
                %/closestDistScale by expected Brownian displacement
                ratioDist2Std = repmat(nnDistTracksE(iTrack)/closestDistScale,...
                    1,timeWindow) ./ dispBrown1;

                %make ratios larger than maxStdMult equal to maxStdMult
                ratioDist2Std(ratioDist2Std > maxStdMult) = maxStdMult;

                %expand end's search radius multiplication factor if possible
                brownStdMultModE = max([brownStdMultModE; ratioDist2Std]);

            end

            %determine the long vectors of the search ellipses at track
            %starts for all time gaps
            longVecS1 = repmat((brownStdMultModS .* dispBrown1 * sqrtDim),probDim,1) .* ...
                repmat(directionMotion,1,timeWindow);

            %determine the long vectors of the search ellipses at track
            %ends for all time gaps
            longVecE1 = repmat((brownStdMultModE .* dispBrown1 * sqrtDim),probDim,1) .* ...
                repmat(directionMotion,1,timeWindow);

            %determine the short vectors at track starts
            shortVecS1 = repmat((brownStdMultModS .* dispBrown1 * sqrtDim),probDim,1) .* ...
                repmat(perpendicular,1,timeWindow);

            %determine the short vectors at track ends
            shortVecE1 = repmat((brownStdMultModE .* dispBrown1 * sqrtDim),probDim,1) .* ...
                repmat(perpendicular,1,timeWindow);

            %get magnitude and direction of both vectors at track starts
            vecMag = sqrt((diag(longVecS1' * longVecS1))'); %magnitude of both vectors
            longVecDir = longVecS1 ./ repmat(vecMag,probDim,1);   %direction of long vector
            shortVecDir = shortVecS1 ./ repmat(vecMag,probDim,1); %direction of short vector
            
            %make sure that magnitude is within allowed range
            vecMag = max([vecMag;repmat(minSearchRadius(kinTypeStart+1),1,timeWindow)]); %compare to minimum
            vecMag = min([vecMag;maxSearchRadius(kinTypeStart+1,:)]); %compare to maximum
            
            %re-calculate both vectors based on modified magnitudes            
            longVecS1 = repmat(vecMag,probDim,1) .* longVecDir; %new long vector
            shortVecS1 = repmat(vecMag,probDim,1) .* shortVecDir; %new short vector

            %construct additional short vectors for 3D problems and save
            %them
            if probDim == 3
                shortVecS13D = repmat(vecMag,probDim,1) .* repmat(perpendicular3D,1,timeWindow);
                shortVecS3D(:,:,iTrack) = shortVecS13D;
            end

            %get magnitude and direction of both vectors at track ends
            vecMag = sqrt((diag(longVecE1' * longVecE1))');  %magnitude of both vectors
            longVecDir = longVecE1 ./ repmat(vecMag,probDim,1);   %direction of long vector
            shortVecDir = shortVecE1 ./ repmat(vecMag,probDim,1); %direction of short vector
            
            %make sure that magnitude is larger than minimum allowed and
            %smaller than maximum allowed
            vecMag = max([vecMag;repmat(minSearchRadius(kinTypeEnd+1),1,timeWindow)]); %compare to minimum
            vecMag = min([vecMag;maxSearchRadius(kinTypeEnd+1,:)]); %compare to maximum
            
            %re-calculate both vectors based on modified magnitudes            
            longVecE1 = repmat(vecMag,probDim,1) .* longVecDir; %new long vector
            shortVecE1 = repmat(vecMag,probDim,1) .* shortVecDir; %new short vector

            %save values for this track
            %             dispBrown(:,iTrack) = dispBrown1;
            longVecS(:,:,iTrack) = longVecS1;
            longVecE(:,:,iTrack) = longVecE1;
            shortVecS(:,:,iTrack) = shortVecS1;
            shortVecE(:,:,iTrack) = shortVecE1;

            %construct additional short vectors for 3D problems and save
            %them
            if probDim == 3
                shortVecE13D = repmat(vecMag,probDim,1) .* repmat(perpendicular3D,1,timeWindow);
                shortVecE3D(:,:,iTrack) = shortVecE13D;
            end

        otherwise

            %take direction of motion to be along x and construct
            %perpendicular(s)
            if probDim == 2
                directionMotion = [1 0]';
                perpendicular = [0 1]';
            else
                directionMotion = [1 0 0]';
                perpendicular = [0 1 0]';
                perpendicular3D = [0 0 1]';
            end

            %calculate the expected displacement along x (= along y) due to
            %brownian motion for all time gaps
            dispBrown1 = undetBrownStd * timeScalingBrown;

            %copy brownStdMult into vector that might be modified using
            %local density
            brownStdMultModS = brownStdMult'; %for track start
            brownStdMultModE = brownStdMult'; %for track end

            %if local density information is used to expand search radius ...
            if useLocalDensity

                %divide the track's nearest neighbor distance/closestDistScale
                %by expected Brownian displacement
                ratioDist2Std = repmat(nnDistTracksS(iTrack)/closestDistScale,...
                    1,timeWindow) ./ dispBrown1;

                %make ratios larger than maxStdMult equal to maxStdMult
                ratioDist2Std(ratioDist2Std > maxStdMult) = maxStdMult;

                %expand start's search radius multiplication factor if possible
                brownStdMultModS = max([brownStdMultModS; ratioDist2Std]);

                %divide the track's nearest neighbor distance/closestDistScale
                %by expected Brownian displacement
                ratioDist2Std = repmat(nnDistTracksE(iTrack)/closestDistScale,...
                    1,timeWindow) ./ dispBrown1;

                %make ratios larger than maxStdMult equal to maxStdMult
                ratioDist2Std(ratioDist2Std > maxStdMult) = maxStdMult;

                %expand end's search radius multiplication factor if possible
                brownStdMultModE = max([brownStdMultModE; ratioDist2Std]);

            end
            
            %determine the long vector of the search ellipse at track
            %starts for all time gaps
            longVecS1 = repmat((brownStdMultModS .* dispBrown1 * sqrtDim),probDim,1) .* ...
                repmat(directionMotion,1,timeWindow);

            %determine the long vector of the search ellipse at track
            %ends for all time gaps
            longVecE1 = repmat((brownStdMultModE .* dispBrown1 * sqrtDim),probDim,1) .* ...
                repmat(directionMotion,1,timeWindow);

            %determine the short vector at track starts
            shortVecS1 = repmat((brownStdMultModS .* dispBrown1 * sqrtDim),probDim,1) .* ...
                repmat(perpendicular,1,timeWindow);

            %determine the short vector at track ends
            shortVecE1 = repmat((brownStdMultModE .* dispBrown1 * sqrtDim),probDim,1) .* ...
                repmat(perpendicular,1,timeWindow);

            %get magnitude and direction of both vectors at track starts
            vecMag = sqrt((diag(longVecS1' * longVecS1))'); %magnitude of both vectors
            longVecDir = longVecS1 ./ repmat(vecMag,probDim,1);   %direction of long vector
            shortVecDir = shortVecS1 ./ repmat(vecMag,probDim,1); %direction of short vector
            
            %make sure that magnitude is larger than minimum allowed and
            %smaller than maximum allowed
            vecMag = max([vecMag;repmat(minSearchRadius(kinTypeStart+1),1,timeWindow)]); %compare to minimum
            vecMag = min([vecMag;maxSearchRadius(kinTypeStart+1,:)]); %compare to maximum
            
            %re-calculate both vectors based on modified magnitudes            
            longVecS1 = repmat(vecMag,probDim,1) .* longVecDir; %new long vector
            shortVecS1 = repmat(vecMag,probDim,1) .* shortVecDir; %new short vector

            %construct additional short vectors for 3D problems
            if probDim == 3
                shortVecS13D = repmat(vecMag,probDim,1) .* repmat(perpendicular3D,1,timeWindow);
                shortVecS3D(:,:,iTrack) = shortVecS13D;
            end
            
            %get magnitude and direction of both vectors at track ends
            vecMag = sqrt((diag(longVecE1' * longVecE1))'); %magnitude of both vectors
            longVecDir = longVecE1 ./ repmat(vecMag,probDim,1);   %direction of long vector
            shortVecDir = shortVecE1 ./ repmat(vecMag,probDim,1); %direction of short vector
            
            %make sure that magnitude is larger than minimum allowed and
            %smaller than maximum allowed
            vecMag = max([vecMag;repmat(minSearchRadius(kinTypeEnd+1),1,timeWindow)]); %compare to minimum
            vecMag = min([vecMag;maxSearchRadius(kinTypeEnd+1,:)]); %compare to maximum
            
            %re-calculate both vectors based on modified magnitudes
            longVecE1 = repmat(vecMag,probDim,1) .* longVecDir; %new long vector
            shortVecE1 = repmat(vecMag,probDim,1) .* shortVecDir; %new short vector

            %save values for this track
            %             dispBrown(:,iTrack) = dispBrown1;
            longVecS(:,:,iTrack) = longVecS1;
            longVecE(:,:,iTrack) = longVecE1;
            shortVecS(:,:,iTrack) = shortVecS1;
            shortVecE(:,:,iTrack) = shortVecE1;

            %construct additional short vectors for 3D problems
            if probDim == 3
                shortVecE13D = repmat(vecMag,probDim,1) .* repmat(perpendicular3D,1,timeWindow);
                shortVecE3D(:,:,iTrack) = shortVecE13D;
            end
            
    end %(switch trackType)

end %(for iTrack = 1 : numTracks)

%%%%% ~~ the end ~~ %%%%%

