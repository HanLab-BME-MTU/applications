function dataList = trajectoryAnalysisMainSeedStates(dataList,distance,time,timePoints,constants)
%does the first fit & classification of trajectory
%

groupLength = constants.MINLENGTH;
PROB2SIDES = constants.PROB2SIDES;
MAXDELETED = constants.MAXDELETED;

if groupLength == 2
    
    
    %------ FIND DELETED TIMEPOINTS
    tpDiff = diff(timePoints);
    %write state: -X where X frames habe been deleted, 0 elsewhere
    state = 1-tpDiff;
    
    %get goodIdx: where less than max have been deleted
    goodIdx = find(state >= -MAXDELETED);
    
    %---END FIND DELETED TIMEPOINTS
    
    %------ FIND SIGNIFICANT DISTANCE CHANGE
    %test: is distance change significantly larger than deltaDSigma?
    %as the variances for the measured values are a priori unknown, we should
    %use a t-test. The number of degrees of freedom, however, is very large, as
    %all pixels of the image count -> use gauss
    
    testValue = abs(dataList(goodIdx,9)./dataList(goodIdx,10)); %deltaD/deltaDSigma (has to be positive for the test)
    testLimit = norminv(PROB2SIDES,0,1); %we use a two-sided test
    
    %H0: delta == 0 / H1: delta ~=0
    testOutcome = testValue>testLimit; %==1 if H1
    
    %-----END FIND SIGNIFICANT DISTANCE CHANGE
    
    
    %------ CLASSIFY DATA
    %growth = 1
    %shrinkage = 2
    %pause = 3
    %undetermined = 0
    %frame deleted = -1, -2 etc, depending on how many frames have been
    %deleted in between two timepoints
    
    %if testOutcome==1&deltaD>0: growth. <: shrinkage. ==0 undetermined (we
    %can't tell about pauses yet). After that, write gaps
    state(goodIdx) = testOutcome.*(1.5-0.5*sign(dataList(goodIdx,9)));
    
    
    %write state into dataList
    dataList(:,3) = state;
    
    %-----END CLASSIFY DATA
    
    %---------CALCULATE DATA
    toCalcIdx = find((dataList(:,4) == 0) & (dataList(:,3) == 1) | (dataList(:,3) == 2));
    toCalcStart = dataList(toCalcIdx,1); %read tp
    toCalcEnd   = dataList(toCalcIdx,2);
    
    %calc deltaD and sigma (sqrt(sum(sigma^2)))
    deltaD = distance(toCalcEnd,1)-distance(toCalcStart,1);
    deltaDSigma = sqrt(distance(toCalcEnd,2).^2 + distance(toCalcStart,2).^2);
    
    time = dataList(toCalcIdx,7);
    
    dataList(toCalcIdx,[4,5,9,10]) = [deltaD./time,deltaDSigma./time,deltaD,deltaDSigma];
    %-----END CALCULATE DATA
    
    
else
    %do triplet etc. fitting here
    disp('this option has not been implemented yet!')
end