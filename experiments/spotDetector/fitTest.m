function [numDist,ncordList,ampList,bg,statistics,debugData]=fitTest(data,cordList,idxList,dataSize,dataProperties,estNoise,DEBUG)
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
F_TEST_PROB=dataProperties.F_TEST_PROB;

MAX_POS_DELTA= [7 7 7];

%init debug parameters
if nargin < 7 | isempty(DEBUG)
    DEBUG = 0;
else
    debugData = struct('exitflag',[],'output',[]);
end

%init output (could remain empty!)
nCordList = [];
ampList = [];
bg = [];
statistics = [];


%shift coords to center
cent=ceil(dataSize/2);
%cordList=cordList-ones(size(cordList,1),1)*cent;

% optim options
numDist = 0;
options = optimset('Jacobian','on','Display','off','TolFun',1e-6,'TolX',1e-6);

% find minimal box for gaussian fit
[ig jg kg]=ind2sub(dataSize,idxList);

%dermine minimal surrounding box
minV=[min(ig) min(jg) min(kg)];
maxV=[max(ig) max(jg) max(kg)];
mskDataSize=maxV-minV+1;

%make odd size
mo=~rem(mskDataSize,2);
mskDataSize=mskDataSize+mo;
gIdxList=sub2ind(mskDataSize, ig-minV(1)+1, jg-minV(2)+1, kg-minV(3)+1);

%shift coord
shiftC=([minV(1)-1, minV(2)-1, minV(3)-1]+ceil(mskDataSize/2));
cordList=cordList-ones(size(cordList,1),1)*shiftC;

%-------------------------------------- FIT N ------------------------------------
% FIT N (N=number of cluster spots found in loc max)
nsp=size(cordList,1);
% free params:
degFree=4*nsp+1;

%setup parms and boundaries
lb=[  0];                         %lb of background
%chJ ub=[255];                      %ub        "
ub = [1]; %-> intensities are transformed onto [0 100]

%get maxData for transform
maxData = max(data(:));

%TO PREVENT BADLY CONDITIONED MATRIX: transform data onto interval [~0 1]
TRANSFACT = 1/maxData;
transData = data*TRANSFACT;

bg=min(transData(:));

sca=(max(transData(:))-min(transData(:)));
parms=[bg];

%position       intensity scaling     rest
for i=1:nsp
    lb=[cordList(i,:)-MAX_POS_DELTA           0   lb];
    ub=[cordList(i,:)+MAX_POS_DELTA 50000  ub];
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
    chi1= sum(Res1(:).^2)/(length(Res1(:))-degFree);
    
    QAll=(gaussgrad'*gaussgrad)^-1;
    
    %=======--TEST FOR SIGNIFICANCE--======
    %for all spots (single or multi) that are fitted: test amplitude whether
    %above noise or not
    %if several spots: test whether distance between them is significant (in
    %pixels). in case this test fails, throw out worst spot, fit with N-1, then
    %still try mixture model
    
    %testDistanceandAmplitudes will return new parms, Q.
    [parms,QAll,deletedSpotNumber,rmIdx] = testDistanceAndAmplitudes(parms,QAll,chi1,dataProperties,0);
    
    if ~isempty(deletedSpotNumber)
        nsp = nsp - length(deletedSpotNumber);
        % update number of degrees of freedom
        degFree=4*nsp+1;
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
    nCt=1;
    for i = 1:nsp
        %boundary
        intensIdx=i*4;
        
        %transform amplitudes again (parms only, because lb and ub have not been transformed back!)
        parms(ampAndBGIdx) =  parms(ampAndBGIdx)*TRANSFACT;
        % for lower bound intensity: use estimated image noise*2
        nlb=[ncordList(i,:) - MAX_POS_DELTA 2*estNoise lb];
        nub=[ncordList(i,:)+MAX_POS_DELTA 50000 ub];
        nlb(intensIdx+4*nCt)=0.2*parms(intensIdx);
        % --- shouldn't we try 2-3 times with random init?
        nparms=[ncordList(i,:)+(2-4*rand(1,3)) 0.5*parms(intensIdx) parms];
        nlb(intensIdx+4*nCt)=0.5*parms(intensIdx);
        
        [nparms,resnorm,dummy,exitflag,output] = lsqnonlin(@distTestError,nparms,nlb,nub,options,transData,gIdxList,mskDataSize,dataProperties);
        
        %transform back parms
        nlp = length(nparms);
        nAmpAndBGIdx = [4:4:nlp,nlp];
        nparms(nAmpAndBGIdx) = nparms(nAmpAndBGIdx)/TRANSFACT;
        
        [nGaussFit, nGaussGrad] = multiGaussFit(mskDataSize,nparms,dataProperties);
        
        ndegFree=4*(nsp+nCt)+1;
        
        %calc chi and Q-matrix
        Res2= nGaussFit(gIdxList)-data;
        chi2= sum(Res2(:).^2)/(length(Res2(:))-ndegFree);
        nQAll=(nGaussGrad'*nGaussGrad)^-1;
        
        if DEBUG
            debugData(end+1).exitflag = exitflag;
            debugData(end).output = output;
            debugData(end).resnorm = resnorm;
            debugData(end).gaussResnorm = sum(Res2(:).^2);
        end
        
        %test the fit if significantly improved
        fValue=(chi1/degFree)/(chi2/ndegFree);
        prob=fcdf(fValue,degFree,ndegFree);
        
        if (prob>F_TEST_PROB)
            %test again whether the spots are significant
            [nparms,nQAll,deletedSpotNumber] = testDistanceAndAmplitudes(nparms,nQAll,chi2,dataProperties,1);
            %if we had to delete anything this time, we do not accept the N+1-fit
            
            if isempty(deletedSpotNumber)
                QAll = nQAll;
                lb=nlb;
                ub=nub;
                parms=nparms;
                ampAndBGIdx = nAmpAndBGIdx;
                nCt=nCt+1;
            else %transform back parms
                parms(ampAndBGIdx) = parms(ampAndBGIdx)/TRANSFACT;
            end
            
        else %transform back parms
            parms(ampAndBGIdx) = parms(ampAndBGIdx)/TRANSFACT;
        end
    end %for-loop
    
    
    %read out parameters - only if nsp>0
    
    %index to postion values
    posIdx=sort([1:4:(4*(nsp+nCt-1)) 2:4:(4*(nsp+nCt-1)) 3:4:(4*(nsp+nCt-1))]);
    ampIdx=[4:4:(4*(nsp+nCt-1))];
    
    %read Q, coordinates and amplitudes
    Q=QAll(posIdx,posIdx);
    qAmp = QAll(ampIdx,ampIdx);
    ncordList=reshape(parms(posIdx),3,nsp+nCt-1)';
    ampList=parms(ampIdx);
    
    
    if nCt>1
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
    numDist=nsp+nCt-1;
    %shift coords back
    ncordList=ncordList+ones(size(ncordList,1),1)*shiftC;
end

%we can always assign the background...
bg=parms(end);