function [numDist,ncordList,ampList,bg,statistics,debugData]=fitTest(data,cordList,idxList,dataSize,dataProperties,DEBUG)
%FITTEST 
%
%SYNOPSIS : [numDist,ncordList,statistics]=fitTest(data,cordList,idxList,dataSize)
%
% INPUT   data       :  vector with raw data
%              cordList :  list with center coordinates of detected spots
%              idxList    : vector of mask
%              dataSize:  original size of data (3D shape)
%
% OUTPUT numDist   : number of underlying gaussians found
%                ncordList : estimated coords
%                ampList   : amplitudes of spots
%                bg           : common background value
%                statistics  : detection statistics (struct: chi squared, snr, multi/single spot,)
%
% Recent changes: - amplitude is tested aginst the background, and joint
% fits are tested for whether they are significantly far form each other (otherwise, we will fit with one spot less)

% c: 7/6/01 dT

%CONST DEFINITIONS
F_TEST_PROB = dataProperties.F_TEST_PROB;

MAX_POS_DELTA = dataProperties.FILTERPRM(4:6);

%init debug parameters
if nargin < 6 || isempty(DEBUG)
    DEBUG = 0;
    debugData = [];
else
    debugData = struct('exitflag',[],'output',[]);
end

%init output (could remain empty!)
ncordList = [];
ampList = [];
bg = [];
statistics = [];


% optim options. We want to terminate the algorithm when it is within 0.001
% pixel.
numDist = 0;
options = optimset('Jacobian','on','Display','off','TolFun',1e-20,'TolX',1e-3);

% find minimal box for gaussian fit
[ig jg kg]=ind2sub(dataSize,idxList);

%dermine minimal surrounding box
minV=[min(ig) min(jg) min(kg)];
maxV=[max(ig) max(jg) max(kg)];
mskDataSize=maxV-minV+1;

% %make odd size. This is some relic from earlier times, when everything
% had to be odd size to work properly. At least, this time the image is not
% made smaller, so I don't mess with it.
mo=~rem(mskDataSize,2);
mskDataSize=mskDataSize+mo;
gIdxList=sub2ind(mskDataSize, ig-minV(1)+1, jg-minV(2)+1, kg-minV(3)+1);


% shift coord. ShiftC is the coordinate of the lower left corner of the
% localImage plus the distance from the lower left corner of the localImage
% to the center of the localImage
localShift = ceil(mskDataSize/2);
shiftC=([minV(1)-1, minV(2)-1, minV(3)-1] + localShift);
cordList=cordList-ones(size(cordList,1),1)*shiftC;

%-------------------------------------- FIT N ------------------------------------
% FIT N (N=number of cluster spots found in loc max)
nsp=size(cordList,1);
% free params:
numFreeParms=4*nsp+1;

%setup parms and boundaries
lb=[  0];                         %lb of background
%chJ ub=[255];                      %ub        "
ub = [1]; %-> intensities are transformed onto [0 1]

%get maxData for transform
maxData = max(data(:));

%TO PREVENT BADLY CONDITIONED MATRIX: transform data onto interval [~0 1]
TRANSFACT = 1/maxData;
transData = data*TRANSFACT;


bg=min(transData(:));

sca=(max(transData(:))-bg);
parms=[bg];

%position       intensity scaling     rest
% lower bound for intensity: 0. Using the estimate noise seemed to be a
% good idea until I found that if lsqnonlin hits a boundary, it still gives
% some kind of reasonable fit - with the lowerbound as amplitude! Of
% course, the artificially high amplitude will pass the statistical tests.
% upper bound for intensity: 1.1 (b/c data is on 0/1)
for i=1:nsp
    lb=[cordList(i,:)-MAX_POS_DELTA, 0,   lb];
    ub=[cordList(i,:)+MAX_POS_DELTA, 1.1,  ub];
    parms=[cordList(i,:) sca parms];
end;

redoN = 1;
%loop until all spots have passed the test
while redoN
    redoN = 0;
    %calculate the fit, return exitflag, output for debugging
    [parms,resnorm,dummy,exitflag,output] = lsqnonlin(@distTestError,parms,lb,ub,options,transData,gIdxList,mskDataSize,dataProperties);
    
    %backtransform data. don't forget background!
    lp = length(parms);
    ampAndBGIdx = [4:4:lp,lp];
    parms(ampAndBGIdx) = parms(ampAndBGIdx)/TRANSFACT;
    
    %calculate the residuals and the gradient
    [gaussFit, gaussgrad] = multiGaussFit(mskDataSize,parms,dataProperties);
    
    Res1=gaussFit(gIdxList)-data;
    degreesOfFreedom = (length(Res1(:))-numFreeParms);
    chi1= sum(Res1(:).^2)/degreesOfFreedom;
    
    QAll=(gaussgrad'*gaussgrad)^-1;
    
    %=======--TEST FOR SIGNIFICANCE--======
    %for all spots (single or multi) that are fitted: test amplitude whether
    %above noise or not
    %if several spots: test whether distance between them is significant (in
    %pixels). in case this test fails, throw out worst spot, fit with N-1, then
    %still try mixture model
    
    %testDistanceandAmplitudes will return new parms, Q.
    [parms,QAll,deletedSpotNumber,rmIdx,...
        debugData(end+1).testValue] = testDistanceAndAmplitudes(...
        parms,QAll,chi1,dataProperties,0,degreesOfFreedom);
    
    if ~isempty(deletedSpotNumber)
        nsp = nsp - length(deletedSpotNumber);
        % update number of free parameters. Degrees of Freedom will be
        % updated after the fitting
        numFreeParms=4*nsp+1;
        %send the data through another fitting loop
        if nsp>0
            redoN = 1;
            parms([4:4:4*nsp,4*nsp+1]) = parms([4:4:4*nsp,4*nsp+1])*TRANSFACT; %prepare parms
            ub(rmIdx) = [];
            lb(rmIdx) = [];
        end
    end
    
    if DEBUG
        debugData(end+1).exitflag = exitflag;
        debugData(end).output = output;
        debugData(end).resnorm = resnorm;
        debugData(end).gaussResnorm = sum(Res1(:).^2);
    end
    
end %while redoN

%index to postion values
posIdx=sort([1:4:(4*nsp) 2:4:(4*nsp) 3:4:(4*nsp)]);
ampIdx=[4:4:(4*(nsp))];

%read new Qcord, ncordList
Q=QAll(posIdx,posIdx);
qAmp = QAll(ampIdx,ampIdx);
ncordList=reshape(parms(posIdx),3,nsp)';

%sprintf('%05.3f \n',parms)


if nsp>0 %do N+1-fit only if there are any spots left!
    %--------------------------------------FIT N+1-------------------------------------
    
    % nCt counts the number of new tags
    nCt=0;
    
    % local image is the image of the jointly fitted tags. Put NaNs
    % wherever there is no data
    localImage = repmat(NaN,mskDataSize);
    localImage(gIdxList) = data;
    [subImageGridY,subImageGridX,subImageGridZ] = ndgrid(...
        [-MAX_POS_DELTA(1):MAX_POS_DELTA(1)],...
     [-MAX_POS_DELTA(2):MAX_POS_DELTA(2)],...
     [-MAX_POS_DELTA(3):MAX_POS_DELTA(3)]);
    
    % loop through all the jointly fitted spots
    for i = 1:nsp
        % loop for every spot while new spots are being found
        failedTest = 0;
        newIdx = [];
        while ~failedTest
        
        
        %boundaries and starting values.
        
        intensIdx=(i+nCt)*4;
        %transform amplitudes again (parms only, because lb and ub have not been transformed back!)
        parms(ampAndBGIdx) =  parms(ampAndBGIdx)*TRANSFACT;
        % for lower bound intensity: be very low: we don't want to
        % introduce artefacts
        nlb=[ncordList(i,:) - MAX_POS_DELTA 0 lb];
        nub=[ncordList(i,:)+MAX_POS_DELTA 2 ub];
        % as first guess for the new parameter we take the centroid of the
        % subimage around the spot. If we are in the second iteration, the
        % patch center will be in between the old and the new spot. This of
        % course assumes that tags will always be arranged in some kind of
        % a triangle, and never in a line with the brightest in the center,
        % but for these cases there's the tracker. 
        
        % for subimage: interpolate around patchCenter. This is slower than
        % cutting out an image via stamp3D, but easier to understand
        patchCenter = mean(ncordList([i+nCt,newIdx],:),1) + localShift;
        % interpolate. Put NaN wherever we get out of range
        subImage = interp3(localImage,subImageGridX + patchCenter(2),...
            subImageGridY + patchCenter(1),...
            subImageGridZ + patchCenter(3),...
            '*linear',NaN);
        centroid = centroid3D(subImage.^10);
        % careful: centroid works in imageCoords - this code is written in
        % matrix-coords, though.
        centroid = centroid([2,1,3]);
        % (centroid - (maxPosDelta + 1)) is the shift relative to the
        % subImage-center. Since especially in the case of kinetochores
        % being eaten by SPBs the centroid is "attracted" by the brighter
        % tag, we go 5x the centerShift.
        centerShift = 5*(centroid - MAX_POS_DELTA - 1);
        % keep 0.5*intensity for the new spot
        nparms = [patchCenter + centerShift - localShift 0.5*parms(intensIdx) parms];
        
        
        [nparms,resnorm,dummy,exitflag,output] = ...
            lsqnonlin(@distTestError,nparms,nlb,nub,...
            options,transData,gIdxList,mskDataSize,dataProperties);
        
        %transform back parms
        nlp = length(nparms);
        nAmpAndBGIdx = [4:4:nlp,nlp];
        nparms(nAmpAndBGIdx) = nparms(nAmpAndBGIdx)/TRANSFACT;
        
        % recalc residuals and gradient
        [nGaussFit, nGaussGrad] = multiGaussFit(mskDataSize,nparms,dataProperties);
        
        % new number of free parameters
        newNumFreeParms=4*(nsp+nCt+1)+1;
        
        %calc chi and Q-matrix
        Res2= nGaussFit(gIdxList)-data;
        newDegreesOfFreedom = (length(Res2(:))-newNumFreeParms);
        chi2= sum(Res2(:).^2)/newDegreesOfFreedom;
        nQAll=(nGaussGrad'*nGaussGrad)^-1;
        
        if DEBUG
            debugData(end+1).exitflag = exitflag;
            debugData(end).output = output;
            debugData(end).resnorm = resnorm;
            debugData(end).gaussResnorm = sum(Res2(:).^2);
        end
        
        %test the fit if significantly improved
        % df???
%         fValue=(chi1)/(chi2);
%         prob=fcdf(fValue,degreesOfFreedom,newDegreesOfFreedom);
        fValue=(chi1)/(chi2);
        prob=fcdf(fValue,numFreeParms,newNumFreeParms);

        % disp(sprintf('%1.4f',prob));
        if (prob>F_TEST_PROB)
            %test again whether the spots are significant
            [nparms,nQAll,deletedSpotNumber,dummy,debugData(end+1).testValue] = ...
                testDistanceAndAmplitudes(...
                nparms,nQAll,chi2,dataProperties,1,newDegreesOfFreedom);
            %if we had to delete anything this time, we do not accept the N+1-fit
            
            if isempty(deletedSpotNumber)
                QAll = nQAll;
                lb=nlb;
                ub=nub;
                parms=nparms;
                ampAndBGIdx = nAmpAndBGIdx;
                nCt=nCt+1;
                posIdx=sort([1:4:(4*(nsp+nCt)) 2:4:(4*(nsp+nCt)) 3:4:(4*(nsp+nCt))]);
                ncordList=reshape(parms(posIdx),3,(nsp+nCt))';
                newIdx = [newIdx, 0] + 1;
                % make sure we compare to the last good fit, not to the
                % initial one
                numFreeParms = newNumFreeParms;
                chi1 = chi2;
            else %transform back parms and quit loop
                failedTest = 1;
                parms(ampAndBGIdx) = parms(ampAndBGIdx)/TRANSFACT;
            end
            
        else %transform back parms
            failedTest = 1;
            parms(ampAndBGIdx) = parms(ampAndBGIdx)/TRANSFACT;
        end
        end % while-loop
    end %for-loop
    
    
    %read out parameters - only if nsp>0
    
    %index to postion values
    posIdx=sort([1:4:(4*(nsp+nCt)) 2:4:(4*(nsp+nCt)) 3:4:(4*(nsp+nCt))]);
    ampIdx=[4:4:(4*(nsp+nCt))];
    
    %read Q, coordinates and amplitudes
    Q=QAll(posIdx,posIdx);
    qAmp = QAll(ampIdx,ampIdx);
    ncordList=reshape(parms(posIdx),3,nsp+nCt)';
    ampList=parms(ampIdx);
    
    
    if nCt>0
        statistics.parms=parms;
        statistics.multi=1;
        statistics.chi=chi2;
    else
        statistics.parms=parms;
        statistics.multi=0;
        statistics.chi=chi1;
    end;
    
    %comp SNR
    statistics.snr=parms(ampIdx)/sqrt(statistics.chi);
    
    
    statistics.Qxx=Q;
    statistics.qAmp = qAmp;
    numDist=nsp+nCt;
    %shift coords back
    ncordList=ncordList+ones(size(ncordList,1),1)*shiftC;
end

%we can always assign the background...
bg=parms(end);
