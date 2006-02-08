function [testRatios]=fitTestFirstRound(data,cordList,idxList,dataSize,dataProperties,DEBUG)
%FITTESTFIRSTROUND fits the data to find the cutoff for the amplitude test.
% It's a combination of fitTest and testDistanceAndAmplitudes

MAX_POS_DELTA = dataProperties.FILTERPRM(4:6);

% %init debug parameters
% if nargin < 6 || isempty(DEBUG)
%     DEBUG = 0;
%     debugData = [];
% else
%     debugData = struct('exitflag',[],'output',[]);
% end

% %init output (could remain empty!)
% ncordList = [];
% ampList = [];
% bg = [];
% statistics = [];


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


%calculate the fit, return exitflag, output for debugging
parms = lsqnonlin(@distTestError,parms,lb,ub,options,transData,gIdxList,mskDataSize,dataProperties);

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

%----------- get amplitude ratio and return
%init parameters
nSpots=floor(length(parms)/4); %number of spots (last parm is bg)
endList = 4*nSpots;

%read amps
ampIdx = [4:4:endList];
amp   = parms( ampIdx );

%read Qa (amps)
Qa = QAll(ampIdx,ampIdx);

% get testRatios
testRatios = amp'./sqrt(diag(Qa)*chi1);
