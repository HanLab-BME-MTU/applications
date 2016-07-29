function [transDiffAnalysisRes,errFlag] = basicTransientDiffusionAnalysis(tracks,...
    probDim,alphaValues,alphaAsym,minDuration,plotRes,confRadMin,checkAsym)
%BASICTRANSIENTDIFFUSIONANALYSIS detects potential diffusion segments of a track and performs MSS analysis on these segments
%
%SYNOPSIS [transDiffAnalysisRes,errFlag] = basicTransientDiffusionAnalysis(tracks,...
%     extractType,probDim,alphaValues,minDuration,plotRes,confRadMin,checkAsym)
%
%INPUT  tracks      : -- EITHER --
%                     Output of trackWithGapClosing (matrix),
%                     -- OR --
%                     Output of trackCloseGapsKalman (structure, possibly
%                     with merges/splits).
%
%       probDim     : Problem dimensionality.
%                     Optional. Default: 2.
%       alphaValues : Alpha-value for classification. Can take the values
%                     0.2, 0.1, 0.05 and 0.01. One can enter one value, in
%                     which case it will be used for both confined and
%                     directed, or two values, where the first will be used
%                     for confined and the second for directed.
%                     Optional. Default: 0.05 for both.
%       minDuration : Minimum category duration. One can enter one value,
%                     in which case it will be used for both confined and
%                     directed, or two values, where the first will be used
%                     for confined and the second for directed.
%                     Optional. Default: [8 2].
%       plotRes     : 1 to plot results, 0 otherwise.
%                     Optional. Default: 0.
%                     Results can be plotted only if problem is 2D.
%                     color-coding:
%                     *blue: confined diffusion.
%                     *cyan: normal diffusion.
%                     *magenta: super diffusion.
%                     *black: unclassified.
%       confRadMin  : 1 to calculate the confinement radius of confined
%                     particles using the minimum positional standard
%                     deviation, 0 to calculate it using the mean
%                     positional standard deviation.
%                     Optional. Default: 0.
%
%OUTPUT transDiffAnalysisRes : And array (size = number of tracks) with the
%                     field ".segmentClass", which contains the fields:
%           .momentScalingSpectrum: (Number of classification
%                     subparts)-by-(20+probDim) matrix, where each row
%                     contains the information for one classification
%                     subpart, and the columns store the following:
%                     (1) Start frame of subpart.
%                     (2) End frame of subpart.
%                     (3) Classification of subpart: 1 = confined, 2 =
%                         free, 3 = directed.
%                     (4) MSS slope resulting in classification.
%                     (5-11) Generalized diffusion coefficients for moment
%                            orders 0-6.
%                     (12-18) Scaling power for moment orders 0-6.
%                     (19) Normal diffusion coefficient (from the MSD).
%                     (20) Confinement radius, if subpart is classified as
%                          confined (NaN otherwise).
%                     (21/22/23) Center of subpart, if subpart is
%                                classified as confined (NaN otherwise).
%           .momentScalingSpectrum1D: same as above but in 1D for asym
%           tracks, additionally includes preferred direction value
%           .asymmetry, .asymmetryAfterMSS: NOT IMPLEMENTED RIGHT NOW.
%
%       errFlag         : 0 if executed normally, 1 otherwise.
%
%REMARKS While tracks do not have to be linear in order to be asymmetric,
%the last analysis step assumes that tracks are linear.
%
%Khuloud Jaqaman, March 2008
%Updated Tony Vega, July 2016

%% Output

transDiffAnalysisRes = [];
errFlag = 0;

%% Input

%check whether tracks were input
if nargin < 1
    disp('--trackTransientDiffusionAnalysis1: Please input at least the tracks to be analyzed!');
    errFlag = 1;
    return
end

if nargin < 2 || isempty(probDim)
    probDim = 2;
end

if nargin < 3 || isempty(alphaValues)
    alphaValues = 0.05;
end

if nargin < 4 || isempty(alphaAsym)
    alphaValues = 0.05;
end

if nargin < 5 || isempty(minDuration)
    minDuration = [8 2];
elseif length(minDuration) == 1
    minDuration = [minDuration minDuration];
end

if nargin < 6 || isempty(plotRes)
    plotRes = 0;
elseif plotRes == 1 && probDim ~= 2
    disp('--trackTransientDiffusionAnalysis1: Cannot plot tracks if problem is not 2D!');
    plotRes = 0;
end

if nargin < 7 || isempty(confRadMin)
    confRadMin = 0;
end

if errFlag
    disp('--trackTransientDiffusionAnalysis1: Please fix input variables');
    return
end


%define window sizes
windowAsym = 5;
windowMSS = 21;
windowMSSMin = 20;
halfWindowAsym = (windowAsym - 1) / 2;
halfWindowMSS = (windowMSS - 1) / 2;

%define duration of a Brownian segment that leads to lumping it with the
%segments around it
windowBrown = 20;

%specify MSS analysis moment orders
momentOrders = 0 : 6;

%% Track extraction for analysis

%store input tracks in a new variable
tracksInput = tracks;

%extract segments for analysis if tracks were input as a structure that
%might contain merges and splits
%the point is to reduce compound tracks that contain merges and splits into
%simple separate tracks
%thus this step is not necessary if the tracks were input as a matrix,
%which by definition does not contain unresolved compound tracks.
if isstruct(tracks)

    %get number of input tracks from structure
    numInputTracks = length(tracksInput);

    clear tracks

    switch extractType

        case 1 %retrieve every track segment separately

            [tracks,dummy,compTrackStartRow,numSegments] = ...
                convStruct2MatIgnoreMS(tracksInput);

        case 2 %make the longest track possible, given all the merges and splits

            disp('Sorry - not implemented yet!')
            errFlag = 1;
            return

    end

else

    %get number of input tracks from matrix
    numInputTracks = size(tracksInput,1);

    %indicate rows where tracks start (trivial in this case)
    compTrackStartRow = (1 : numInputTracks)';

    %indicate number of segments in each track (1 for all tracks)
    numSegments = ones(numInputTracks,1);

end

%get number of track segments to be analyzed
numTrackSegments = size(tracks,1);

%get track segment start times, end times and life times
trackSEL = getTrackSEL(tracks);

%find track segments that are long enough for analysis

%     indx4analysis = find(trackSEL(:,3) >= windowMSSMin);
    indx4analysis = find(trackSEL(:,3) >= windowAsym);


indxNot4analysis = setdiff((1:numTrackSegments)',indx4analysis);

%% Rolling window classification

%reserve memory %This will only have certain traits, probably none of these
trackSegmentClassRes = repmat(struct('asymmetry',NaN(1,3),...
    'momentScalingSpectrum',NaN(1,21+probDim),...
    'momentScalingSpectrum1D',NaN(1,21+probDim),...
    'asymmetryAfterMSS',NaN(1,3)),...
    numTrackSegments,1);

%go over all analyzable track segments
gaussDeriv = cell(numTrackSegments,1); 
j=1;
for iTrack = indx4analysis'

    %% Asymmetry analysis to get directed parts of the track segment

    %get track segment start, end and life times
    trackSELCurrent = trackSEL(iTrack,:);

    %Just initialize for now
    trackClassAsym = [trackSELCurrent(1:2) NaN];


    %% Maximum displacement analysis to divide track segment into parts with potentially different diffusion behavior

    %assign initial MSS classification
    %all parts classified above as asymmetric get a -1
    trackClassMSS = trackClassAsym;
    trackClassMSS(:,3) = -trackClassMSS(:,3);
    trackClassMSS(trackClassMSS(:,3)~=-1,3) = NaN;
    trackClassMSS(:,4:20+probDim) = NaN;
    oldPart2newPartMap = (1 : size(trackClassMSS,1))';

    %find the length of each part
    trackPartLength = trackClassMSS(:,2) - trackClassMSS(:,1) + 1;

    %find track segment parts that are not classified as asymmetric and
    %that are longer than windowMSSMin
    trackParts2analyze = find(trackPartLength >= windowMSSMin);
    %if there are parts to analyze ...
    if ~isempty(trackParts2analyze)

        %go over these parts ...
        for iPart = trackParts2analyze'

            %get starting point of this part
            trackPartStart = trackClassMSS(oldPart2newPartMap(iPart),1);

            %get number of MSS analysis rolling windows in this part
            numRollWindows = trackPartLength(iPart) - windowMSS + 1;

            %if number of rolling windows is larger than the minimum
            %required duration of a classification, proceed with rolling
            %window analysis
            if numRollWindows > minDuration(1)

                %initialize max displacement vector
                maxDisplacement = NaN(numRollWindows,1); 
                %go over all windows
                for iWindow = 1 : numRollWindows

                    %get window start and end points
                    startPoint = trackPartStart + iWindow - 1;
                    endPoint = startPoint + windowMSS - 1;
                    
                    test = tracks(iTrack,8*(startPoint-1)+1:8*endPoint);
                    xTest = test(1:8:end);
                    yTest = test(2:8:end);
                    X =[xTest',yTest'];
                    D = pdist(X,'euclidean');
                    maxDisplacement(iWindow) = max(D);
                    

                    
                end %(for iWindow = 1 : numRollWindows)
                %% Derivative of Maximum Displacement
                h =1;
                y = maxDisplacement;
                [out] = filterGauss1D(y, 2, 'symmetric'); 
                der = diff(out)/h; 
                der = abs(der); 
                normDer = der./(maxDisplacement(1:end-1)+maxDisplacement(2:end));
                gaussDeriv{iTrack} = normDer;

            else
                gaussDeriv{iTrack} = [];
            end
        end
    end
    j=j+1;
end
clear trackClassMSS


for iTrack = indx4analysis'
    %% Separate tracks that change and those that don't
    trackFull = gaussDeriv{iTrack};
    if ~isempty(trackFull)
        thresh = max(trackFull); %level values derived from simulation
        if thresh >=0.15
            level = 0.05;
        else
            level= 0.02;
        end
 
        %% Asymmetry detection first
        %1. Check track to see if any sections have asymmetry
        %this classification scheme is taken from Huet et al (BJ 2006)
        %it classifies tracks as asymmetric or not, based on the scatter of
        %positions along them
        
        if checkAsym       
            %Divide track into segments using asym minimum
            [segPointsA] = findDiffSegments(trackFull,level,halfWindowMSS,windowAsym);% using halfWindowMSS to look at same part of track analyzed above
            n = 1:length(segPointsA)-1;
            difference = segPointsA(n+1)-segPointsA(n);
            trackSELCurrent = trackSEL(iTrack,:);
            partClassAsym = NaN(length(n),3);
            %go over all of these segments and determine if they are
            %asymmetric
            for k = 1:length(segPointsA)-1 
            startPoint = trackSELCurrent(1) +segPointsA(k)-1;
            endPoint  = startPoint+ difference(k)-1;

                %get the particle positions along the track
                coordX = (tracks(iTrack,8*(startPoint-1)+1:8:8*endPoint))';
                coordY = (tracks(iTrack,8*(startPoint-1)+2:8:8*endPoint))';
                %coordZ = (tracks(iTrack,8*(startPoint-1)+3:8:8*endPoint))';%Not implemented yet
                coordXY = [coordX coordY];

                %determine whether the track is sufficiently asymmetric
                [~,asymFlag] = asymDeterm2D3D(coordXY(:,1:probDim),alphaAsym);

                %classify track as ...
                %1 = linear, if the asymmetry parameter is larger than the threshold
                %0 = not linear, if the asymmetry parameter is smaller than the
                %threshold
                %otherwise, keep track classification as undetermined
                partClassAsym(k,1) = startPoint;
                partClassAsym(k,2) = endPoint;
                partClassAsym(k,3) = asymFlag;
            end


        %find indices of all tracks classified as asymmetric
        indxAsym = find(partClassAsym(:,3) == 1);
        else
            %Other values which are not defined, given value of NaN or
            %similar
                indxAsym = 0;
                [segPointsA] = findDiffSegments(trackFull,level,halfWindowMSS,windowMSSMin);% halfWindowMSS,windowMSSMin
        end
    
        %% Back to normal
        %Initialize all variables that will be saved later
        n = 1:length(segPointsA)-1;
        difference = segPointsA(n+1)-segPointsA(n);
        pointClassMSS = NaN(length(n),1);
        mssSlope = pointClassMSS;
        normDiffCoef = pointClassMSS;
        confRadTmp = pointClassMSS;
        centerTmp = NaN(length(n),probDim);
        genDiffCoef = NaN(length(n),length(momentOrders));
        scalingPower = NaN(length(n),length(momentOrders));
        partClassMSS = NaN(length(n),3);
        
        partClassMSS1D = NaN(length(n),3);
        mssSlope1D = pointClassMSS;
        genDiffCoef1D = NaN(length(n),length(momentOrders));
        scalingPower1D = NaN(length(n),length(momentOrders));
        normDiffCoef1D = pointClassMSS;
        confRadius1D = NaN(length(n),2);
        prefDir = NaN(length(n),2);
        trackCenter = NaN(length(n),2);
        %% Diffusion Analysis
        %Now go through segments and get diffusion classification.
        trackSELCurrent = trackSEL(iTrack,:);
        for k = 1:length(segPointsA)-1 
        startPoint = trackSELCurrent(1) +segPointsA(k)-1;
        endPoint  = startPoint+ difference(k)-1;
        
            if ismember(k,indxAsym)    

                        %get the positions in this track and their standard deviations

                    [pointClass,mssSlopeT,genDiffCoefT,scalingPowerT,normDiffCoefT,trackCenter(k,:),confRadius1D(k,:),prefDir(k,:)] = asymmetricDiffusion(startPoint,endPoint,alphaAsym,probDim,tracks,iTrack);

                    %since not all track segments are linear, put analysis results in their
                    %proper place among all track segment
                    partClassMSS1D(k,1)= startPoint;
                    partClassMSS1D(k,2)= endPoint;              
                    partClassMSS1D(k,3) = partClassAsym(k,3);
                    partClassMSS(k,1)= startPoint;
                    partClassMSS(k,2)= endPoint;
                    partClassMSS(k,3)= pointClass;
                    mssSlope1D(k) = mssSlopeT;
                    genDiffCoef1D(k,:) = genDiffCoefT;
                    scalingPower1D(k,:) = scalingPowerT;
                    normDiffCoef1D(k) = normDiffCoefT;

            else
                        [pointClassMSS(k),mssSlope(k),genDiffCoef(k,:),scalingPower(k,:),normDiffCoef(k)] = trackMSSAnalysis(...
                            tracks(iTrack,8*(startPoint-1)+1:8*endPoint),...
                            probDim,momentOrders,alphaValues(1));


                        if pointClassMSS(k) <= 1
                        [confRadTmp(k),centerTmp(k,:)] = estimConfRad(tracks(iTrack,8*(startPoint-1)+1:8*endPoint),probDim,confRadMin);
%{
                                    if length(startPoint:endPoint) <= 100 %Could even go up to 100, look into further
                                        xCoord = (tracks(iTrack,8*(startPoint-1)+1:8:8*endPoint))';
                                        yCoord = (tracks(iTrack,8*(startPoint-1)+2:8:8*endPoint))';
                                        [tClassMSS]  = immobileDetection(xCoord,yCoord,mssSlope(k));
                                        %% here
                                        pointClassMSS(k) = tClassMSS;
                                    end
    %}
                        else
                            confRadTmp(k) = NaN;
                            centerTmp(k,:) = NaN(1,probDim);    
                        end 
                        partClassMSS1D(k,1)= startPoint;
                        partClassMSS1D(k,2)= endPoint;
                        partClassMSS(k,1)= startPoint;
                        partClassMSS(k,2)= endPoint;
                        partClassMSS(k,3)= pointClassMSS(k);
            end
                    
        end
        
         %% Reclassify if segments next to each other are the same
                            %merge subparts that now have the same classification
                            partClassTmp = partClassMSS;
                            pointClassMSSAsymT = partClassMSS1D(:,3);
                    for iSubpart = 1 : length(pointClassMSS)-1
                        iSubpartPlus1 = iSubpart + 1;
                        while ( (iSubpartPlus1 <= length(pointClassMSS)) && ...
                                ( (pointClassMSS(iSubpart) == pointClassMSS(iSubpartPlus1)) || ...
                                (isnan(pointClassMSS(iSubpart)) && isnan(pointClassMSS(iSubpartPlus1))) ) && ...
                                ( (pointClassMSSAsymT(iSubpart) == pointClassMSSAsymT(iSubpartPlus1)) || ... %same asym
                                (isnan(pointClassMSSAsymT(iSubpart)) && isnan(pointClassMSSAsymT(iSubpartPlus1))) ) )
                            
                            partClassMSS1D(iSubpart,2) = partClassMSS1D(iSubpartPlus1,2);
                            partClassMSS1D(iSubpartPlus1,1) = partClassMSS1D(iSubpart,1);
                            
                            partClassMSS(iSubpart,2) = partClassMSS(iSubpartPlus1,2);
                            partClassMSS(iSubpartPlus1,1) = partClassMSS(iSubpart,1);
                            iSubpartPlus1 = iSubpartPlus1 + 1;
                        end
                    end
                    [~,uniqueParts] = unique(partClassMSS(:,1));
                    partClassMSS = partClassMSS(uniqueParts,:);
                    partClassMSS1D = partClassMSS1D(uniqueParts,:);
                    
                    %Also check if there are symmetric unclassified segments
                    unclassCheck  = isnan(partClassMSS(:,3)).*isnan(partClassMSS1D(:,3));
                    if sum(unclassCheck)>0 && size(partClassMSS,1)>1
                        [partClassMSS,partClassMSS1D] = mergeUnclassSegments(unclassCheck,partClassMSS,partClassMSS1D,trackFull,halfWindowMSS);
                    end
                    
                    %if parts have been merged, redo diffusion analysis on
                    %these parts. Otherwise continue
                    if size(partClassMSS,1) < size(partClassTmp,1)
                        numSeg = size(partClassMSS,1);
                        pointClassMSS = NaN(numSeg,1);
                        partClassAsym = NaN(numSeg,1);
                        mssSlope = pointClassMSS;
                        normDiffCoef = pointClassMSS;
                        confRadTmp = pointClassMSS;
                        centerTmp = NaN(numSeg,probDim);
                        genDiffCoef = NaN(numSeg,length(momentOrders));
                        scalingPower = NaN(numSeg,length(momentOrders));
                        
                                mssSlope1D = pointClassMSS;
                        genDiffCoef1D = NaN(numSeg,length(momentOrders));
                        scalingPower1D = NaN(numSeg,length(momentOrders));
                        normDiffCoef1D = pointClassMSS;
                        confRadius1D = NaN(numSeg,2);
                        prefDir = NaN(numSeg,2);
                        trackCenter = NaN(numSeg,2);
                        
                        for k = 1:numSeg 
                        startPoint = partClassMSS(k,1);
                        endPoint  = partClassMSS(k,2);
                        alphaAsym = 0.05;%alphaValues(2);

                                    %get the particle positions along the track
                                    coordX = (tracks(iTrack,8*(startPoint-1)+1:8:8*endPoint))';
                                    coordY = (tracks(iTrack,8*(startPoint-1)+2:8:8*endPoint))';
%                                     coordZ = (tracks(iTrack,8*(startPoint-1)+3:8:8*endPoint))';
                                    coordXY = [coordX coordY];

                                    %determine whether the track is sufficiently asymmetric
                                    
                                        [~,asymFlag] = asymDeterm2D3D(coordXY(:,1:probDim),alphaAsym);

                                    %classify track as ...
                                    %1 = linear, if the asymmetry parameter is larger than the threshold
                                    %0 = not linear, if the asymmetry parameter is smaller than the
                                    %threshold
                                    %otherwise, keep track classification as undetermined

                                    if asymFlag ==1
                                        
                                        [pointClass,mssSlopeT,genDiffCoefT,scalingPowerT,normDiffCoefT,trackCenter(k,:),confRadius1D(k,:),prefDir(k,:)] = asymmetricDiffusion(startPoint,endPoint,alphaAsym,probDim,tracks,iTrack);

                                        %since not all track segments are linear, put analysis results in their
                                        %proper place among all track segment
                                        partClassMSS1D(k,1)= startPoint;
                                        partClassMSS1D(k,2)= endPoint;              
                                        partClassMSS1D(k,3) = asymFlag;
                                        partClassMSS(k,1)= startPoint;
                                        partClassMSS(k,2)= endPoint;
                                        partClassMSS(k,3)= pointClass;

                                        mssSlope1D(k) = mssSlopeT;
                                        genDiffCoef1D(k,:) = genDiffCoefT;
                                        scalingPower1D(k,:) = scalingPowerT;
                                        normDiffCoef1D(k) = normDiffCoefT;
                                    else

                                        [pointClassMSS(k),mssSlope(k),genDiffCoef(k,:),scalingPower(k,:),normDiffCoef(k)] = trackMSSAnalysis(...
                                            tracks(iTrack,8*(startPoint-1)+1:8*endPoint),...
                                            probDim,momentOrders,alphaValues(1));
     
                                        if pointClassMSS(k) <= 1  
                                        [confRadTmp(k),centerTmp(k,:)] = estimConfRad(tracks(iTrack,8*(startPoint-1)+1:8*endPoint),probDim,confRadMin);
                                        %{
                                            if length(startPoint:endPoint) <= 100 %Could even go up to 100, look into further
                                                xCoord = (tracks(iTrack,8*(startPoint-1)+1:8:8*endPoint))';
                                                yCoord = (tracks(iTrack,8*(startPoint-1)+2:8:8*endPoint))';
            %                                     center  = centerTmp(k,:);
                                                [tClassMSS]  = immobileDetection(xCoord,yCoord,mssSlope(k));
                                                %% Here
                                                pointClassMSS(k) = tClassMSS;
                                            end
                                        %}
                                        else
                                            confRadTmp(k) = NaN;
                                            centerTmp(k,:) = NaN(1,probDim);
                                        end
                                        partClassMSS(k,3) = pointClassMSS(k);
                                   end
                            
                            
                        end
                    
                    end

                  
    else
        %analyze whole track 
        trackSELCurrent = trackSEL(iTrack,:);
        
        %Fill in empty 1D data
        mssSlope1D = NaN(1,1);
        mssSlope = NaN(1,1);
        genDiffCoef1D = NaN(1,length(momentOrders));
        genDiffCoef = NaN(1,length(momentOrders));
        scalingPower1D = NaN(1,length(momentOrders));
        scalingPower = NaN(1,length(momentOrders));
        normDiffCoef1D = NaN(1,1);
        normDiffCoef = NaN(1,1);
        confRadius1D = NaN(1,2);
        confRadTmp = NaN;
        centerTmp = NaN(1,probDim);
        
        prefDir = NaN(1,2);
        trackCenter = NaN(1,2);
        partClassMSS1D = [trackSELCurrent(:,1:2) NaN];
        
        if checkAsym ==1

            %get the particle positions along the track
            coordX = (tracks(iTrack,1:8:end))';
            coordY = (tracks(iTrack,2:8:end))';
%             coordZ = tracks(iTrack,3:8:end)';
            coordXY = [coordX coordY];

            %determine whether the track is sufficiently asymmetric
            [~,asymFlag] = asymDeterm2D3D(coordXY(:,1:probDim),alphaAsym);

            if asymFlag ==1

                [pointClass,mssSlope1D,genDiffCoef1D,scalingPower1D,normDiffCoef1D,trackCenter,confRadius1D,prefDir] = asymmetricDiffusion(trackSELCurrent(:,1),trackSELCurrent(:,2),alphaAsym,probDim,tracks,iTrack);

                partClassMSS1D = [trackSELCurrent(:,1:2) asymFlag];
                partClassMSS = [trackSELCurrent(:,1:2) pointClass];

            else
                
                trackSELCurrent = trackSEL(iTrack,:);
                [pointClassMSS,mssSlope,genDiffCoef,scalingPower,normDiffCoef] = trackMSSAnalysis(...
                    tracks(iTrack,:),...
                    probDim,momentOrders,alphaValues(1));
                %Get confinement info if necessary
                if pointClassMSS <= 1
                    [confRadTmp,centerTmp] = estimConfRad(tracks(iTrack,:),probDim,confRadMin);
                else
                    confRadTmp = NaN;
                    centerTmp = NaN(1,probDim);
                end 
                partClassMSS =[trackSELCurrent(:,1:2) pointClassMSS]; 
                partClassMSS1D = [trackSELCurrent(:,1:2) asymFlag];

            end
        else
                [pointClassMSS,mssSlope,genDiffCoef,scalingPower,normDiffCoef] = trackMSSAnalysis(...
                    tracks(iTrack,:),...
                    probDim,momentOrders,alphaValues(1));
        %Get confinement info if necessary
                if pointClassMSS <= 1
                    [confRadTmp,centerTmp] = estimConfRad(tracks(iTrack,:),probDim,confRadMin);
                else
                    confRadTmp = NaN;
                    centerTmp = NaN(1,probDim);
                end 
                partClassMSS =[trackSELCurrent(:,1:2) pointClassMSS]; 
                    
        end
                trackClassAsym = partClassMSS;
                trackClassAsym(:,3) =NaN;

    end
    
    %% Store somewhere
    trackClassMSS = [partClassMSS,mssSlope,genDiffCoef,scalingPower, normDiffCoef,confRadTmp,centerTmp];
    trackClassMSS1D = [partClassMSS1D,mssSlope1D,genDiffCoef1D,scalingPower1D, normDiffCoef1D,confRadius1D,trackCenter,prefDir]; %COMPLETE
    trackSegmentClassRes(iTrack).asymmetry = partClassMSS1D;% [1 352 NaN]
    trackSegmentClassRes(iTrack).momentScalingSpectrum = trackClassMSS;
    trackSegmentClassRes(iTrack).momentScalingSpectrum1D = trackClassMSS1D;
    trackSegmentClassRes(iTrack).asymmetryAfterMSS = trackClassAsym;% [1 352 NaN] Verify whether this is needed
    
end
%% Store trivial nonclassification information for tracks that are not classifiable
if ~isempty(indxNot4analysis')
for iTrack = indxNot4analysis' 
    trackSELCurrent = trackSEL(iTrack,:);
    trackSegmentClassRes(iTrack).asymmetry(1:2) = trackSELCurrent(1:2);
    trackSegmentClassRes(iTrack).momentScalingSpectrum(1:2) = trackSELCurrent(1:2);
     trackSegmentClassRes(iTrack).momentScalingSpectrum1D(1:2) = trackSELCurrent(1:2);
    trackSegmentClassRes(iTrack).asymmetryAfterMSS(1:2) = trackSELCurrent(1:2);
end
end
%% save results in output structure

%reserve memory
segmentClass = struct('asymmetry',[],'momentScalingSpectrum',[],'momentScalingSpectrum1D',[],'asymmetryAfterMSS',[]);
transDiffAnalysisRes = repmat(struct('segmentClass',segmentClass),numInputTracks,1);

%go over all input tracks
for iTrack = 1 : numInputTracks

    %go over the segments of each track
    for iSegment = 1 : numSegments(iTrack)

        %store the segment's classification results
        transDiffAnalysisRes(iTrack).segmentClass(iSegment,1).asymmetry = ...
            trackSegmentClassRes(compTrackStartRow(iTrack)+iSegment-1).asymmetry;
        transDiffAnalysisRes(iTrack).segmentClass(iSegment,1).momentScalingSpectrum = ...
            trackSegmentClassRes(compTrackStartRow(iTrack)+iSegment-1).momentScalingSpectrum;
        transDiffAnalysisRes(iTrack).segmentClass(iSegment,1).momentScalingSpectrum1D = ...
            trackSegmentClassRes(compTrackStartRow(iTrack)+iSegment-1).momentScalingSpectrum1D;
        transDiffAnalysisRes(iTrack).segmentClass(iSegment,1).asymmetryAfterMSS = ...
            trackSegmentClassRes(compTrackStartRow(iTrack)+iSegment-1).asymmetryAfterMSS;

    end %(for iSegment = 1 : numSegments(iTrack))

end %(for iTrack = 1 : numInputTracks)

%% plotting

%plot results if requested
if plotRes
    plotTracksTransDiffAnalysis2D(tracksInput,transDiffAnalysisRes,[],1,[],[],checkAsym,[]);
end



end

%% Subfunctions
function  [partClassMSS,partClassMSS1D] = mergeUnclassSegments(unclassCheck,partClassMSS,partClassMSS1D,trackFull,halfWindowMSS)
    %cycle through segments
    check = find(unclassCheck);
%         [~,check] =sort(difference(difference < windowMSSMin)); 
        while ~isempty(check)
        %If short segment is first, connect to succeeding
            if check(1) == 1
                partClassMSS(1,2) = partClassMSS(2,2);
                partClassMSS1D(1,2) = partClassMSS1D(2,2);
                partClassMSS(1,3) = partClassMSS(2,3);
                partClassMSS1D(1,3) = partClassMSS1D(2,3);
                partClassMSS(2,:) =[];
                partClassMSS1D(2,:) =[];
                unclassCheck  = isnan(partClassMSS(:,3)).*isnan(partClassMSS1D(:,3));
                check = find(unclassCheck);
                continue
            end
             %If segment is last, connect to preceding
            if check(1) == length(unclassCheck) %changed from end
                partClassMSS(end-1,2) = partClassMSS(end,2);
                partClassMSS1D(end-1,2) = partClassMSS1D(end,2);
                partClassMSS(end,:) =[];
                partClassMSS1D(end,:) =[];
                unclassCheck  = isnan(partClassMSS(:,3)).*isnan(partClassMSS1D(:,3));
                check = find(unclassCheck);
                continue
            end

            %If segment is somewhere in middle
            if ~isempty(check)
                %If preceeding and succeeding segments are the same
                %diffusion, merge everything
                if partClassMSS(check(1)+1,3) == partClassMSS(check(1)-1,3) && isnan(partClassMSS1D(check(1)+1,3)) == isnan(partClassMSS1D(check(1)-1,3))
                    partClassMSS(check(1)-1,2) = partClassMSS(check(1)+1,2);
                    partClassMSS1D(check(1)-1,2) = partClassMSS1D(check(1)+1,2);
                    partClassMSS(check(1)+1,:) =[];
                    partClassMSS1D(check(1)+1,:) =[];
                    partClassMSS(check(1),:) =[];
                    partClassMSS1D(check(1),:) =[];
                elseif partClassMSS1D(check(1)+1,3) == partClassMSS1D(check(1)-1,3) && isnan(partClassMSS(check(1)+1,3)) == isnan(partClassMSS(check(1)-1,3))
                    partClassMSS(check(1)-1,2) = partClassMSS(check(1)+1,2);
                    partClassMSS1D(check(1)-1,2) = partClassMSS1D(check(1)+1,2);
                    partClassMSS(check(1)+1,:) =[];
                    partClassMSS1D(check(1)+1,:) =[];
                    partClassMSS(check(1),:) =[];
                    partClassMSS1D(check(1),:) =[];  
                 %Otherwise score the segments and merge to weaker score
                else
                    if trackFull(partClassMSS(check(1)+1,1)-(partClassMSS(1,1)+halfWindowMSS)) > trackFull(partClassMSS(check(1),1)-(partClassMSS(1,1)+halfWindowMSS))
                        partClassMSS(check(1)-1,2) = partClassMSS(check(1),2);
                        partClassMSS1D(check(1)-1,2) = partClassMSS1D(check(1),2);
                        partClassMSS(check(1),:) =[];
                        partClassMSS1D(check(1),:) =[];
                    else
                        partClassMSS(check(1),2) = partClassMSS(check(1)+1,2);
                        partClassMSS1D(check(1),2) = partClassMSS1D(check(1)+1,2);
                        partClassMSS(check(1)+1,:) =[];
                        partClassMSS1D(check(1)+1,:) =[];
                    end
                end
                       
                unclassCheck  = isnan(partClassMSS(:,3)).*isnan(partClassMSS1D(:,3));
                check = find(unclassCheck);
            end
        end
end
function [pointClass, mssSlopeT,genDiffCoefT,scalingPowerT,normDiffCoefT,trackCenter,confRadius1D,prefDir] = asymmetricDiffusion(startPoint,endPoint,alphaAsym,probDim,tracks, iTrack)
                    trackCoordX = tracks(iTrack,8*(startPoint-1)+1:8:8*endPoint)';
                    deltaCoordX = tracks(iTrack,8*(startPoint-1)+5:8:8*endPoint)';
                    trackCoordY = tracks(iTrack,8*(startPoint-1)+2:8:8*endPoint)';
                    deltaCoordY = tracks(iTrack,8*(startPoint-1)+6:8:8*endPoint)';
                    trackCoordZ = tracks(iTrack,8*(startPoint-1)+3:8:8*endPoint)';
                    deltaCoordZ = tracks(iTrack,8*(startPoint-1)+7:8:8*endPoint)';
                    trackCoord = [trackCoordX trackCoordY trackCoordZ];
                    deltaCoord = [deltaCoordX deltaCoordY deltaCoordZ];
                    trackCoord = trackCoord(:,1:probDim);
                    deltaCoord = deltaCoord(:,1:probDim);

                    %project onto direction of motion
                    [posAlongDir,deltaPosAlongDir] = projectCoordOntoDir(trackCoord,...
                        deltaCoord,[],[]);

                    %construct matrix of linear tracks with projected positions
                    trackCoord2 = [posAlongDir zeros(length(posAlongDir),3) deltaPosAlongDir zeros(length(posAlongDir),3)]';
                    trackCoord2 = trackCoord2(:)';
            %         iAsym = iAsym + 1;
            %         tracksAsym(iAsym,:) = trackCoord2;
                    momentOrders = 0 : 6;
                    [pointClass,mssSlopeT,genDiffCoefT,scalingPowerT,normDiffCoefT] = ...
                    trackMSSAnalysis(trackCoord2,1,momentOrders,alphaAsym);
                %% Confinement of asym 
                    xyCoord = [trackCoordX trackCoordY];

                    %find the eignevalues of the variance-covariance matrix of this track's
                    %positions
                    [eigenVec,eigenVal] = eig(nancov(xyCoord(:,1:probDim)));
                    eigenVal = diag(eigenVal);

                    %calculate the confinement radius along the preferred direction of
                    %motion
                    confRadius1D(1,2) = sqrt( max(eigenVal) * (3) );

                    %calculate the confinement radius perpendicular to the preferred
                    %direction of motion
                    confRadius1D(1,1) = sqrt( mean(eigenVal(eigenVal~=max(eigenVal))) * (probDim + 1) );

                    %calculate the track's center
                    trackCenter = nanmean(xyCoord(:,1:probDim));

                    %store the preferred direction of motion
                    prefDir = eigenVec(:,eigenVal==max(eigenVal))';
end

function [confRadTmp,centerTmp] = estimConfRad(tracks,probDim,confRadMin)

%get subpart's coordinates
xCoord = tracks(1:8:end)';
yCoord = tracks(2:8:end)';
zCoord = tracks(3:8:end)';
xyzCoord = [xCoord yCoord zCoord];

%find the eigenvalues and eigenvectors of the variance-covariance
%matrix of this track's positions
eigenVal = eig(nancov(xyzCoord(:,1:probDim)));

%calculate the track's confinement radius
if confRadMin
    confRadTmp = sqrt( min(eigenVal) * (probDim + 2) );
else
    confRadTmp = sqrt( mean(eigenVal) * (probDim + 2) );
end

%calculate the track's center
centerTmp = nanmean(xyzCoord(:,1:probDim));
end
function [tClassMSS]  = immobileDetection(xCoord,yCoord,mssSlope)
        %% New part calculate the track's center
        
%         principalComponents = [-0.5536 08334; 0.8334 0.5526];
        principalComponents = [-0.2923  0.9561;  0.9561 0.2923];
        traj = [xCoord,yCoord];
        trajLength = length(traj);
        
            trajMean = nanmean(traj(:,1:2));
            distXYConf = traj- repmat(trajMean,trajLength,1);
            distConf = sqrt(sum(distXYConf.^2,2));
            posXYConf = distXYConf > 0;
            numSwitchConf = length(find(sum(diff(posXYConf),2) > 0));
            fracSwitchConf = numSwitchConf / (trajLength-1);
            sTest(:,1) = fracSwitchConf;
            sTest(:,2)= mssSlope;
            pc = sTest*principalComponents;
            if pc(1) <0%prev 0
                tClassMSS = 0;
            else
                tClassMSS = 1;
            end
        
end
function [segPoints] = findDiffSegments(trackFull,level,halfWindowMSS,windowMSSMin)
        [peaks,locs] = findpeaks(trackFull);%If this slow, try find(trackFull(n) > I(n+1) & trackFull(n) > I(n-1));
%         idxSwitch = find(peaks >= level); %change to level prctile(trackFull,90)/ level
        
        %Get segments and verify that each is a minimum length, if not then
        %join with adjacent segment with smallest peak. May need to repeat
        %after checking that segments are long enough. Make sure entire
        %track is retained 
        peakThresh = peaks(peaks >= level);
        segPoints = [1;locs(peaks >= level);length(trackFull)];
        segPoints(2:end)=segPoints(2:end)+1+halfWindowMSS; %half of windowMin, 1 for deriv offset
        segPoints(end)=segPoints(end)+1+halfWindowMSS; %half of windowMin, 1 for edn correction
        n = 1:length(segPoints)-1;
        difference = segPoints(n+1)-segPoints(n);
% %         difference(1) = difference(1) +1; %offset correction, verify
        checkT = find(difference < windowMSSMin);
% %     if length(checkT)< length(difference)
        ind = difference(difference < windowMSSMin);
        check = checkT(ind  == min(ind));
%         [~,check] =sort(difference(difference < windowMSSMin)); 
        while ~isempty(check)
        %If short segment is first, connect to succeeding
            if check(1) == 1
                segPoints(2,:) =[];
                peakThresh(1) =[];
                
                n = 1:length(segPoints)-1;
                difference = segPoints(n+1)-segPoints(n);
% %                 difference(1) = difference(1) +1; %offset correction, verify
                checkT = find(difference < windowMSSMin);
                ind = difference(difference < windowMSSMin);
                check = checkT(ind  == min(ind));
%                 [~,check] =sort(difference(difference < windowMSSMin));
                continue
            end
             %If segment is last, connect to preceding
            if check(1) == length(difference) %changed from end
                segPoints(end-1,:) = [];
                peakThresh(end) = [];
                                                
                n = 1:length(segPoints)-1;
                difference = segPoints(n+1)-segPoints(n);
% %                 difference(1) = difference(1) +1; %offset correction, verify
                checkT = find(difference < windowMSSMin);
                ind = difference(difference < windowMSSMin);
                check = checkT(ind  == min(ind));
%                 [~,check] =sort(difference(difference < windowMSSMin));
                continue
            end

            %If segment is somewhere in middle, connect to segment with lowest
            %peak
            if ~isempty(check)
                if peakThresh(check(1)) > peakThresh(check(1)-1)
                    segPoints(check(1),:) =[]; %if higher peak is in front, delete peak behind
                    peakThresh(check(1)-1)=[];
                else
                    segPoints(check(1)+1,:) =[];
                    peakThresh(check(1)) =[];
                end
                       
                n = 1:length(segPoints)-1;
                difference = segPoints(n+1)-segPoints(n);
                checkT = find(difference < windowMSSMin);
                ind = difference(difference < windowMSSMin);
                check = checkT(ind  == min(ind));
%                 [~,check] =sort(difference(difference < windowMSSMin));
            end
        end
end
                    %{
                    %% Lower threshold Diffusion Analysis
                    %If any track has been classified as confined or
                    %immobile, check again with lower threshold
                    recheck = ismember(partClassMSS(:,3),[1,0]);
                    if sum(recheck) == 0
%                         continue
                    else
                        
                        tmpPartClassMSS = partClassMSS(recheck,:);
                        for m = 1:size(tmpPartClassMSS,1)
                          testPoints = segPoints2(segPoints2 >= tmpPartClassMSS(m,1) & segPoints2 <= tmpPartClassMSS(m,2)+1);
                          n = 1:length(testPoints)-1;
                          difference = testPoints(n+1)-testPoints(n);
                          pointClassMSS = NaN(length(n),1);
                          mssSlopeTmp = pointClassMSS;
                          normDiffCoefTmp = pointClassMSS;
                          confRadLow = pointClassMSS;
                          centerLow = NaN(length(n),probDim);
                          genDiffCoefTmp = NaN(length(n),length(momentOrders));
                          scalingPowerTmp = NaN(length(n),length(momentOrders));
                          partClassMSSTmp = NaN(length(n),3);
                          trackSELCurrent = trackSEL(iTrack,:);
                                  
                            for k = 1:length(testPoints)-1 
                            startPoint = trackSELCurrent(1) +testPoints(k)-1;
                            endPoint  = startPoint+ difference(k)-1;

                                        [pointClassMSS(k),mssSlopeTmp(k),genDiffCoefTmp(k,:),scalingPowerTmp(k,:),normDiffCoefTmp(k)] = trackMSSAnalysis(...
                                            tracks(iTrack,8*(startPoint-1)+1:8*endPoint),...
                                            probDim,momentOrders,alphaValues);

                                        if pointClassMSS(k) <= 1
                                            [confRadLow(k),centerLow(k,:)] = estimConfRad(tracks(iTrack,8*(startPoint-1)+1:8*endPoint),probDim,confRadMin);
                                            xCoord = (tracks(iTrack,8*(startPoint-1)+1:8:8*endPoint))';
                                            yCoord = (tracks(iTrack,8*(startPoint-1)+2:8:8*endPoint))';
                                            center  = centerLow(k,:);
% %                                             [pointClassMSS(k)]  = immobileDetection(xCoord,yCoord,center);
                                        else
                                            confRadLow(k) = NaN;
                                            centerLow(k,:) = NaN(1,probDim);
                                        end 

                                        partClassMSSTmp(k,1)= startPoint;
                                        partClassMSSTmp(k,2)= endPoint;
                                        partClassMSSTmp(k,3)= pointClassMSS(k);
                                        %Rest of MSS information or not... jsut fill in track

                            end
                            % Connect like segments again
                            %partClassTmp = partClassMSSTmp;
                                for iSubpart = 1 : length(pointClassMSS)-1
                                    iSubpartPlus1 = iSubpart + 1;
                                    
                                    while ( (iSubpartPlus1 <= length(pointClassMSS)) && ...
                                    ( (pointClassMSS(iSubpart) == pointClassMSS(iSubpartPlus1)) || ...
                                    (isnan(pointClassMSS(iSubpart)) && isnan(pointClassMSS(iSubpartPlus1))) ) )
                                    partClassMSSTmp(iSubpart,2) = partClassMSSTmp(iSubpartPlus1,2);
                                    partClassMSSTmp(iSubpartPlus1,1) = partClassMSSTmp(iSubpart,1);
                                    iSubpartPlus1 = iSubpartPlus1 + 1;
                                    end
                                end
                                
                            [~,uniqueParts] = unique(partClassMSSTmp(:,1));
                            partClassMSSTmp = partClassMSSTmp(uniqueParts,:);
                            
                            %create no switch option, nothing changed
                            
                            if sum(ismember(partClassMSSTmp(:,3),2))>0
                                partClassMSSTest =partClassMSS;
%                                 continue
                            elseif ismember(partClassMSSTmp,partClassMSS,'rows')
                                partClassMSSTest =partClassMSS;
%                                 continue
                            else
                                if m ==1
                                    oldClass =partClassMSS(~recheck,:);
%                                     mssSlope0=mssSlope(~recheck,:);
%                                     genDiffCoef0=genDiffCoef(~recheck,:);
%                                     scalingPower0 =scalingPower(~recheck,:);
%                                     normDiffCoef0=normDiffCoef(~recheck,:);
%                                     confRadTmp0=confRadTmp(~recheck,:);
%                                     centerTmp0=centerTmp(~recheck,:);
                                else
                                        oldClass = partClassMSS;

%                                     mssSlope0=mssSlope;
%                                     genDiffCoef0=genDiffCoef;
%                                     scalingPower0 =scalingPower;
%                                     normDiffCoef0=normDiffCoef;
%                                     confRadTmp0=confRadTmp;
%                                     centerTmp0=centerTmp;
                                end
%                                 clear partClassMSS
                                partClassMSSTest = [oldClass;partClassMSSTmp];
%                                 mssSlope=[mssSlope0;mssSlopeTmp];
%                                 genDiffCoef=[genDiffCoef0;genDiffCoefTmp];
%                                 scalingPower =[scalingPower0;scalingPowerTmp];
%                                 normDiffCoef=[normDiffCoef0;normDiffCoefTmp];
%                                 confRadTmp=[confRadTmp0;confRadLow];
%                                 centerTmp=[centerTmp0;centerLow];
                               %Sort out the rest of the variables being
                               %saved
                               %Eventually accomodate for many segments to
                               %be re-analyzed
                            end
                            
                        end
                       partClassMSS = partClassMSSTest; 
                    end
                    
                     %% II Reclassify again if segments next to each other are the same
                     %merge subparts that now have the same classification
                            partClassTmp = partClassMSS;
                            pointClassMSS = partClassMSS(:,3); 
                    for iSubpart = 1 : length(pointClassMSS)-1
                        iSubpartPlus1 = iSubpart + 1;
                        while ( (iSubpartPlus1 <= length(pointClassMSS)) && ...
                                ( (pointClassMSS(iSubpart) == pointClassMSS(iSubpartPlus1)) || ...
                                (isnan(pointClassMSS(iSubpart)) && isnan(pointClassMSS(iSubpartPlus1))) ) )
                            partClassMSS(iSubpart,2) = partClassMSS(iSubpartPlus1,2);
                            partClassMSS(iSubpartPlus1,1) = partClassMSS(iSubpart,1);
                            iSubpartPlus1 = iSubpartPlus1 + 1;
                        end
                    end
                    [~,uniqueParts] = unique(partClassMSS(:,1));
                    partClassMSS = partClassMSS(uniqueParts,:);
                    %Final Diffusion characteristics

                        numSeg = size(partClassMSS,1);
                        pointClassMSS = NaN(numSeg,1);
                        mssSlope = pointClassMSS;
                        normDiffCoef = pointClassMSS;
                        confRadTmp = pointClassMSS;
                        centerTmp = NaN(numSeg,probDim);
                        genDiffCoef = NaN(numSeg,length(momentOrders));
                        scalingPower = NaN(numSeg,length(momentOrders));
                        
                        
                        for k = 1:numSeg 
                        startPoint = partClassMSS(k,1);
                        endPoint  = partClassMSS(k,2);

                                    [pointClassMSS(k),mssSlope(k),genDiffCoef(k,:),scalingPower(k,:),normDiffCoef(k)] = trackMSSAnalysis(...
                                        tracks(iTrack,8*(startPoint-1)+1:8*endPoint),...
                                        probDim,momentOrders,alphaValues);

                                    if pointClassMSS(k) <= 1
                                        [confRadTmp(k),centerTmp(k,:)] = estimConfRad(tracks(iTrack,8*(startPoint-1)+1:8*endPoint),probDim,confRadMin);
                                        xCoord = (tracks(iTrack,8*(startPoint-1)+1:8:8*endPoint))';
                                        yCoord = (tracks(iTrack,8*(startPoint-1)+2:8:8*endPoint))';
                                        center  = centerTmp(k,:);
% %                                         [pointClassMSS(k)]  = immobileDetection(xCoord,yCoord,center);
                                    else
                                        confRadTmp(k) = NaN;
                                        centerTmp(k,:) = NaN(1,probDim);
                                    end 
                            partClassMSS(k,3) = pointClassMSS(k);
                        end
                    
%}