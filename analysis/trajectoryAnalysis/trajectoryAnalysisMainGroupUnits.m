function dataListG = trajectoryAnalysisMainGroupUnits(dataList,distance,time,constants,compatible)
%TRAJECTORYANALYSISMAINGROUPUNITS groups consecutive units of the same type according to the selected strategy 
%trajectoryAnalysisMainGroupUnits(dataListS,distance,time,tp2dataIdx,constants,compatible)
%...update help when done...
%
%



%read constants
PROB2SIDES = constants.PROB2SIDES;
STRATEGY = constants.STRATEGY;
PROBOUTLIER = constants.PROBOUTLIER; %the lower the less accepted (the more outliers)
PROBF = constants.PROBF;%constants.TESTPROBABILITY;

DEBUG = [];%1,2,both

%in the case we stop the evaluation: assign output already
dataListG = dataList;

%get size of dataList
sizeDataList = size(dataList);

%---------FIND GROUPS TO FIT 

%we are interested in all compatible states
states = compatible(:);
%minLength of a group has to be 2 original units
minLength = 2;
%we do not care about how many different
minDifferent = 1;

%find groups
[group2List]=findGroups(dataList(:,3),states,compatible,minLength,minDifferent);

%----END FIND GROUPS TO FIT

%test that we want to fit at all
if isempty(group2List)
    return
end


%--------SWITCH ACCORDING TO STRAGEGY
%STRATEGY 1: start with complete fit, then fit smaller groups as long as not
%            significant or down at original units
%STRATEGY 2: start with two-units fit. Increase length of one unit at a
%            time, until nothing can be enlarged anymore
%

if any(states==3) %testing for pauses, states is [(-X),0,3]
    testPause = 1;
else
    testPause = 0;
end

if DEBUG
    testPause
    collectP = [];
end


switch STRATEGY
    case 1 %fit topDown: try to fit everything, if this does not work, shorten the interval
        
        %------INIT STRATEGY
        
        %init while loop: assign control and storage variables
        nGroups  = size(group2List,1);
        
        ntp2fit = zeros(nGroups,1);
        
        %sort the group2List first
        group2List = sortrows(group2List,[1,2]);
        
        %fill ntp2fit: since cols 1&2 of dataList directly point to the
        %data, we just count the number of indices
        ntp2fit = dataList(group2List(:,2),2)-dataList(group2List(:,1),1)+1;
        
        
        %--END INIT STRATEGY
        
        %go through this loop as long as there are any groups left to fit
        while nGroups > 0
            
            groups2fit = [];
            
            %check whether we can combine groups
            %...put into subfun findAndResolveConflicts...
            validIdx = 1:nGroups;
            checkList = group2List(:,1:2);
            rmGroupNumber = [];
            while ~isempty(validIdx)
                startConflict = validIdx(1);
                validTmp = validIdx(2:end);
                cut = 1;
                %one could conflict with many. look for last
                %conflicting, make this new start until no more conflict
                %if isequal: two groups share a unit (not a timepoint!)
                while any(checkList(startConflict,2)>=checkList(validTmp,1)) %loop will break when validTmp isempty
                    
                    %find which groups directly conflict with current one
                    lastIdx = max(find(checkList(startConflict,2)>=checkList(validTmp,1)));
                    
                    %remove groupsWithinGroups
                    rmGroupNumber = [rmGroupNumber, validTmp(1:lastIdx-1)]; %if lastIdx = 1, nothing is removed
                    
                    %assign new startConflict, cut back validTmp
                    cut = cut+lastIdx;
                    if lastIdx < length(validTmp) %if isequal, all the rest of validTmp is taken
                        startConflict = validTmp(lastIdx);
                        validTmp = validTmp(lastIdx+1:end); %if lastIdx is 1, this makes 2:end
                    else
                        startConflict = [];
                        validTmp = [];
                    end
                    
                end
                
                conflictIdx = validIdx(1:cut);
                
                %update validIdx for the loop to continue
                validIdx = validTmp; %valid-cut
                
                
                %resolve conflicts for groups
                if length(conflictIdx) == 1
                    %this group is completely isolated. Fit.
                    groups2fit = [groups2fit;conflictIdx];
                else
                    %several groups. Join them all in one large group and
                    %fit with the maximum ntp of the groups. Below we make
                    %sure that incompatible units are not fitted at all.
                    
                    %get length of groups
                    conflictingNtp = ntp2fit(conflictIdx);
                    %find maximum number of timepoints
                    ntpMax = max(conflictingNtp);
                    
                    %add a new group at the end of the lists, remove the
                    %merged indices after this loop
                    
                    %remember which to eliminate
                    rmGroupNumber = [rmGroupNumber,conflictIdx];
                    
                    %update ntp
                    ntp2fit(end+1,1) = ntpMax;
                    
                    %update nGroups
                    nGroups = nGroups + 1;
                    
                    %update groups2fit
                    groups2fit = [groups2fit;nGroups];
                    
                    %update group2List
                    group2List(end+1,:) = [group2List(conflictIdx(1),1), group2List(conflictIdx(end),2)];
                    
                end %if length(conflictIdx) == 1
                
            end %while ~isempty(validIdx)
            
            
            %delete merge-entries. Straightforward for all except
            %groups2fit
            if ~isempty(rmGroupNumber)
                rmGroupNumber = rmGroupNumber(:);
                
                tmp = zeros(nGroups,1);
                tmp(groups2fit) = 1;
                tmp(rmGroupNumber) = [];
                groups2fit = find(tmp);
                
                ntp2fit(rmGroupNumber)      = [];
                group2List(rmGroupNumber,:) = [];
                nGroups = nGroups - length(rmGroupNumber);
                
                %since we'll be looking at last/next groups: sort
                [group2List,sortIdx] = sortrows(group2List,[1,2]);
                ntp2fit = ntp2fit(sortIdx);
                %group2List = group2ListOld(sortIdx,:) -> what pointed t
                %g2lo-1, points to sI(1) now
                groups2fit = sortIdx(groups2fit);
            end
            
            
            
            %init new fit-storage
            %newFits: [1:group#,2:startTp,3:endTp,...
            %          4:slope,5:slopeSigma,6:slopeSigmaNull,...
            %          7:chi2Lin,8:nTP,9:meanOrIntercept,10:chi2Mean,11:maxOutlier4Tstats]
            nfListIncrSize = [300, 0];
            nfListSize     = [300, 11];
            newFits        = zeros(nfListSize);
            nfct           = 1;
            
            %loop to do individual jobs, like e.g. fitting.
            %LSCOV is based on c functions while blockdiag to build one
            %huge functional matrix is not.
            for ng = groups2fit';
                
                %calc data that is valid for the whole group
                dataIdx = [dataList(group2List(ng,1),1):dataList(group2List(ng,2),2)]';
                ntp = ntp2fit(ng); %number of (valid) tp in this fit
                shorten = length(dataIdx) - ntp;
                
                for j = 0:shorten
                    
                    %calculate dataIdxFit: which indices to use for this fit
                    dataIdxFit   = dataIdx(1+j:j+ntp); %=1+j:1+j-1+ntp
                    
                    %find the list states of the fit (is there a growth, shrinkage?)
                    
                    %first group in dataList that applies
                    dLFirst = min([find(dataList(:,1)>dataIdxFit(1));sizeDataList(1)+1])-1; %make sure we get no [] here
                    dLLast  = max([find(dataList(:,2)<dataIdxFit(end));0])+1; %make sure we get no [] here
                    
                    %list all states that appear
                    stateList = unique(dataList(dLFirst:dLLast,3));
                    
                    if any(stateList == 1) & any(stateList == 2)
                        %we do not fit at all, counter is not updated
                    else
                        
                        %create functional matries A=dF/dU, B=Y, V=sigmaY^2
                        A = ones(ntp,2); %allocate A
                        A(:,2) = time(dataIdxFit,1); %fill in time
                        B = distance(dataIdxFit,1);
                        Qllii = distance(dataIdxFit,2).^2; %diagonal elements of the covariance matrix
                        V = diag(Qllii); %covariance matrix
                        weightMatrix = diag(1./Qllii); %weight matrix
                        
                        
                        %calculate linear fit. a0+a1*x=y
                        [X,XSigma,XSigmaZeroH] = myLscov(A,B,V);
                        
                        %calculate weighted mean if testpause
                        if testPause
                            %assign weight 1 to the measurement with smallest error
                            wVect = (repmat(min(Qllii,[],1),ntp,1)./Qllii);
                            %calc weightedMean : each dataPoint is multiplied by the corresponding weight, the sum is divided
                            %by the sum of the weights
                            sumWeights = sum(wVect,1);
                            weightedMean = sum(wVect.*B,1)./sumWeights;
                            resMean  = abs(B - repmat(weightedMean,[ntp,1])); %we take abs for the test later
                            chi2Mean = (resMean' * weightMatrix * resMean)/(ntp-1);
                        else
                            chi2Mean = 0; %init this!
                        end
                        
                        %check what new state we would have. If it's not
                        %compatible with an individual state, discard
                        fitState = 1.5-0.5*sign(X(2));
                        if ~testPause && (any(stateList == 1) & fitState == 2) || (any(stateList == 2) & fitState == 1)
                            %do not go further, do not update counter
                        else
                            
                            %calculate chi2
                            res  = abs(A*X-B); %abs here for the outlier testing
                            chi2 = sum(res'*weightMatrix*res)/(ntp-2); %weighted stats!
                            
                            %calculate sigmaApriori
                            sigmaApriori = XSigma/sqrt(XSigmaZeroH);
                            
                            %CHECK FOR OUTLIERS: If there are individual
                            %datapoints that are not compatible with the
                            %fit, discard the fit.
                            %-> because tcdf is SLOW, we just store the
                            %maximum of res./sqrt(Qvvii): if one outlier,
                            %we discard
                            
                            switch testPause
                                case 0 %test linear fit
                                    Qvvii = abs((ones(ntp,1)-A*X)).*Qllii; %diagonal elements of the covariance matrix of the residuals
                                    %pValue = tcdf(res./sqrt(Qvvii),ntp-2);
                                    maxOutlier4T = max(res./sqrt(Qvvii));
                                    meanOrIntercept = A(1,:)*X;
                                    
                                    
                                case 1 %test mean (pause)
                                    %Am = ones(ntp,1), X = weightedMean
                                    Qvvii = abs((1-repmat(weightedMean,1)))*Qllii; %we use the non-normed weights
                                    %pValue = tcdf(resMean./sqrt(Qvvii),ntp-1);
                                    maxOutlier4T = max(res./sqrt(Qvvii));
                                    meanOrIntercept = weightedMean;
                                    
                                    
                            end
                            
                                %fill list - check whether we exceed the list-length
                                if nfct > nfListSize(1)
                                    nfTmp             = newFits; %store old
                                    nfListSize        = nfListSize + nfListIncrSize; %increase size by 1 unit
                                    newFits           = zeros(nfListSize); %init new list
                                    newFits(1:nfct-1,:) = nfTmp; %write back tmp data
                                    clear('nfTmp');
                                end
                                
                                %%newFits: [1:group#,2:startIdx,3:endIdx,...
                                %          4:slope,5:slopeSigma,6:slopeSigmaNull,...
                                %          7:chi2New,8:nTP,9:meanOrIntercept,10:chi2Mean]
                                newFits(nfct,:) = [ng, dataIdxFit(1), dataIdxFit(end),...
                                        X(2), XSigma(2), sigmaApriori(2),...
                                        chi2, ntp, meanOrIntercept, chi2Mean, maxOutlier4T];
                                
                                %update counter
                                nfct = nfct + 1;
                            
                        end   %if ~testPause && (any(stateList == 1) & fitState == 2) && (any(stateList == 2) & fitState == 1)  
                    end %any(stateList == 1) & any(stateList == 2)
                    
                end %for j = 0:shorten :: fill toutDouxList
                
            end %for ng = 1:groups2fit' :: select # of fits, fill newFits
            
            
            
            
            %fill all additional information into list of fits
            
            %remove superfluous entries in newFits & add more needed cols
            newFits(nfct:end,:) = []; %nfct has already been augmented by 1!
            
            
            %set up fits2test to pass to the conflict solver
            %%fits2test: [1:group#,2:startIdx,3:endIdx,...
            %          4:slope,5:slopeSigma,6:slopeSigmaNull,...
            %          7:chi2New,8:nTP,9:meanOrIntercept,10:chi2Mean...
            %          11:pValue,12:maxOutlier4T->pOutlier,...
            %          13: state,14:valid,15:winner]
            fits2resolve = [newFits(:,1:10),zeros(nfct-1,1),newFits(:,11),zeros(nfct-1,3)];
            
            %add p-values, test significance
            %first, test for outliers, then test the rest
            
            if testPause 
                
                %add p-value for outliers
                fits2resolve(:,12) = tcdf(fits2resolve(:,12),fits2resolve(:,8)-1);
                
                %add p-value for test whether the linear fit is significantly better than the mean (Fisher-test)
                fits2resolve(:,11) = fcdf(fits2resolve(:,10)./(fits2resolve(:,7)/2),1,2);
                
                %linear fit is significantly better if f2r>probf => reject
                %pause if it is so & accept only if 
                fits2resolve(:,14) = (fits2resolve(:,11) < PROBF) & (fits2resolve(:,12) < PROBOUTLIER);
                
                if any(DEBUG ==1) & ~isempty(fits2resolve)
                    fits2resolve(:,[2 3 11 12 14])
                end
                
                %write state
                validIdx = find(fits2resolve(:,14));
                fits2resolve(validIdx,13) = 3;
                
            else %test linear fit
                
                %add p-value for outliers
                fits2resolve(:,12) = tcdf(fits2resolve(:,12),fits2resolve(:,8)-2);
                
                %add p-value for test if slope is significantly different
                %from zero 
                fits2resolve(:,11) = tcdf(abs(fits2resolve(:,4)./fits2resolve(:,5)),fits2resolve(:,8)-2);
                
                %accept if pValue > testProbability & no outliers
                fits2resolve(:,14) = (fits2resolve(:,11) > PROB2SIDES) & (fits2resolve(:,12) < PROBOUTLIER);
                
                if any(DEBUG ==2) & ~isempty(fits2resolve)
                    fits2resolve(:,[2 3 11 12 14])
                end
                
                %write state
                validIdx = find(fits2resolve(:,14));
                fits2resolve(validIdx,13) = 1.5-0.5*sign(fits2resolve(validIdx,4));
                
            end
            
            
            %---START CONFLICTS
            
            %resolveConflicts(fits2resolve,groupList)
            %find conflicting fits & resolve
            %...put into subfunction later...
            %start at first valid fit. find any other valid fit starting before
            %endTP1 with the same group. Rinse, repeat.
            
            checkList = fits2resolve(:,2:3);
            while ~isempty(validIdx)
                startConflict = validIdx(1);
                validTmp = validIdx(2:end);
                cut = 1;
                %one could conflict with many. look for last
                %conflicting, make this new start until no more conflict
                %actually, we could just look at group numbers!
                while any(checkList(startConflict,2)>checkList(validTmp,1)) %loop will break when validTmp isempty
                    lastIdx = max(find(checkList(startConflict,2)>checkList(validTmp,1)));
                    %assign new startConflict, cut back validTmp
                    cut = cut+lastIdx;
                    if lastIdx < length(validTmp) %if isequal, all the rest of validTmp is taken
                        startConflict = validTmp(lastIdx);
                        validTmp = validTmp(lastIdx+1:end); %if lastIdx is 1, this makes 2:end
                    else
                        startConflict = [];
                        validTmp = [];
                    end
                end
                
                conflictIdx = validIdx(1:cut);
                
                %update validIdx for the loop to continue
                validIdx = validTmp; %valid-cut
                
                %...subfun resolveConflicts
                
                
                
                %now that we have a conflicting group: resolve the
                %conflict. If the group is of length one, the lone
                %entry is an immediate winner
                
                if length(conflictIdx) == 1
                    %we have a winner
                    fits2resolve(conflictIdx,15) = 1;
                else
                    %we need to look somewhat further
                    
                    %...in other version we have to check for
                    %improvement of fit or p-value if we go from 2 to 3
                    %however, in this STRATEGY, all conflicting fits
                    %have the same length, so we can compare chi2
                    %directly...
                    
                    %loop. since we're doing topdown, indirectly
                    %conflicting fits will be ok, too
                    wct = 1;
                    while  ~isempty(conflictIdx)
                        %the winner is: the one with min.chi2
                        if testPause
                            chiIdx = 10;
                        else
                            chiIdx = 7;
                        end
                        [dummy,winningIdx] = min(fits2resolve(conflictIdx,chiIdx)); 
                        
                        fits2resolve(conflictIdx(winningIdx),15) = wct;
                        %what remains? All that conflict only
                        %indirectly
                        conflictIdxW = conflictIdx(winningIdx);
                        conflictIdx(winningIdx) = [];
                        remainingIdx = find(...
                            fits2resolve(conflictIdxW,2)>=fits2resolve(conflictIdx,3)|fits2resolve(conflictIdxW,3)<=fits2resolve(conflictIdx,2));
                        conflictIdx = conflictIdx(remainingIdx);
                        wct = wct + 1;
                    end
                end
            end
            
            %---END CONFLICTS
            
            %-------ACCEPT WINNERS
            
            %after all the fits, the ntp2fit of the fitted groups decreases
            %by 1
            ntp2fit(groups2fit) = ntp2fit(groups2fit)-1;
            
            %                 
            %write winners into dataListG
            %and update lists
            rmGroupNumber = [];
            allWinnerIdx  = find(fits2resolve(:,15));
            
            for ng = groups2fit' %check only for fitted groups!
                %get index into dataList
                dataListIdx = [group2List(ng,1):group2List(ng,2)]';
                
                %get index into fits2resolve-list
                fitsOfGroup = find(fits2resolve(:,1)==ng);
                
                %find winners of this group. winnerIdx points to fits2resolve                
                winnerIdx = intersect(fitsOfGroup,allWinnerIdx);
                
                if isempty(winnerIdx)
                    %check wheter we can continue fitting is done below
                    %                     if ntp2fit(ng) <= constants.MINLENGTH;
                    %                         rmGroupNumber = [rmGroupNumber;ng];
                    %                     else
                    %                         %don't worry
                    %                     end
                else %we have to update dataListG and all the groupLists (the winning entry will be removed)
                    
                    %remove the list entry
                    rmGroupNumber = [rmGroupNumber;ng];
                    
                    
                    
                    %fill in dataListG first, remember the entries we'll
                    %remove
                    winningEntries = [];
                    
                    for win = 1:length(winnerIdx)
                        %find winning entry in fits2resolve
                        f2rRow = fits2resolve(winnerIdx(win),:);
                        
                        %get first and last entry in dataListG
                        %corresponding to startTP and endTP of this
                        %winner
                        dlgStartIdx = find(dataListG(:,1) == f2rRow(2));
                        dlgEndIdx   = find(dataListG(:,2) == f2rRow(3));
                        
                        dlgSR = dataListG(dlgStartIdx,:);
                        dlgER = dataListG(dlgEndIdx,  :);
                        
                        %build new entry in dataListGroup. 
                        %1:startIdx, 2:endIdx, 3:state, 4:slope, 5:slopeSigma, 6:slopeSigmaAPR,
                        %7:deltaT, 8:(deltaTSigma), 9:deltaD, 10:deltaDSigma, 11:startDistance
                        %don't update deltaT,deltaD - we to this at the
                        %very end only
                        dlgNew = [dlgSR(1),dlgER(2),f2rRow(13),f2rRow(4:6),zeros(1,4),f2rRow(9)];
                        
                        if testPause %make sure there is no growth during a pause
                            dlgNew([4:6])=0; 
                            %startDist, sigma
                            [dlgNew(11),dlgNew(10)] = weightedStats(distance(dlgSR(1):dlgER(2),1),distance(dlgSR(1):dlgER(2),2));
                        end 
                        
                        %remove n-1 entries in dataListGroup
                        dataListG(dlgStartIdx:dlgEndIdx-1,:)=[];
                        %write new entry into last remaining old entry
                        dataListG(dlgStartIdx,:) = dlgNew;
                        
                        %add all the winning dataListEntries to the
                        %winningEntries
                        dlStartIdx = find(dataList(:,1) == f2rRow(2));
                        dlEndIdx   = find(dataList(:,2) == f2rRow(3));
                        winningEntries = [winningEntries;[dlStartIdx:dlEndIdx]']; 
                    end 
                    
                    %update the lists. we look at the remaining entries,
                    %group them and then test for all groups what to do with
                    %them.
                    
                    %remove the winningEntries from the dataListIndices
                    %remainingIdx points to entries in dataListIdx
                    [dummy, remainingIdx] = setdiff(dataListIdx,winningEntries);
                    
                    %next: we group the remainingIdx to find out which
                    %chunks of dataListIdx remain
                    testVector = zeros(max(remainingIdx),1);
                    testVector(remainingIdx)=1;
                    %remainingIdx in testVector form consecutive groups of 1's
                    [remainSE] = findGroups(testVector,1);
                    
                    if isempty(remainSE)
                        %if there were conflicting groups, we have to chop
                        %them now
                        
                        %the previous group can go max to dataListIdx(1)-1,
                        %since these indices point to the dataList that
                        %lists units
                        if ng > 1
                            prevGroupEnd = group2List(ng-1,2);
                            if prevGroupEnd >= dataListIdx(1)
                                group2List(ng-1,2) = dataListIdx(1)-1;
                                try %it could be that the winning group is 1:5 and the conflicting 1:3 -> we would get index 0
                                    ntp2fit(ng-1,1) = min(ntp2fit(ng-1),...
                                        dataList(group2List(ng-1,2),2)-dataList(group2List(ng-1,1),1)+1);
                                catch
                                    ntp2fit(ng-1,1) = 0;
                                end
                            end
                        end
                        %the previous group can start min at dataListIdx(end)+1,
                        %since these indices point to the dataList that
                        %lists units
                        if ng < nGroups
                            nextGroupStart = group2List(ng+1,1);
                            if nextGroupStart <= dataListIdx(end)
                                group2List(ng+1,1) = dataListIdx(end)+1;
                                try %if we go outside dataList here, there is definitely nothing to fit
                                    ntp2fit(ng+1,1) = min(ntp2fit(ng+1),...
                                        dataList(group2List(ng+1,2),2)-dataList(group2List(ng+1,1),1)+1);
                                catch
                                    ntp2fit(ng+1,1) = 0;
                                end
                            end
                        end
                        
                    else
                        
                        
                        %now loop through the the chunks and check for the first
                        %and the last one whether they interfere with another
                        %fit. If yes, we update that fit. If no, and for all
                        %other chunks we test wheter they are longer than
                        %minLength and append them to the lists
                        
                        for chunk = 1:size(remainSE,1)
                            dlStartIdx = dataListIdx(remainSE(chunk,1));
                            dlEndIdx   = dataListIdx(remainSE(chunk,2));
                            
                            %a group is ok if it uses at least two entries in
                            %dataList
                            isValid =  dataList(dlEndIdx) > dataList(dlStartIdx);
                            
                            if chunk == 1 && ng > 1
                                %check whether we interfere with the previous
                                %group
                                prevGroupEnd = group2List(ng-1,2);
                                if prevGroupEnd >= dlStartIdx
                                    %the group ng-1 will end at max at the
                                    %dlEndIdx-1
                                    if prevGroupEnd >= dlEndIdx
                                        %there will be no new chunk, but we have
                                        %to be careful that the prevGroup does
                                        %not get into trouble
                                        isValid = 0;
                                        
                                        group2List(ng-1,2) = dlEndIdx;
                                        
                                        %update list
                                        ntp2fit(ng-1,1) = min(ntp2fit(ng-1),...
                                            dataList(group2List(ng-1,2),2)-dataList(group2List(ng-1,1),1)+1);
                                        
                                        %we do not check for minLength: we have
                                        %to do it at the very end
                                        
                                    else
                                        %there could be a new chunk. prevGroup
                                        %does not change
                                        dlStartIdx = prevGroupEnd+1;
                                        
                                        %check again
                                        isValid = (dlStartIdx < dlEndIdx);
                                        
                                    end
                                else
                                    %test valid below
                                end
                            end
                            
                            if chunk == size(remainSE,1) && ng < nGroups
                                %check interference with next group
                                nextGroupStart = group2List(ng+1,1);
                                if nextGroupStart <= dlEndIdx
                                    %the group ng+1 will start at min at the
                                    %dlStartIdx +1
                                    if nextGroupStart <= dlStartIdx
                                        %there will be no new chunk, but we have
                                        %to be careful that the nextGroup does
                                        %not get into trouble
                                        isValid = 0;
                                        
                                        group2List(ng+1,1) = dlStartIdx;
                                        
                                        %update list
                                        ntp2fit(ng+1,1) = min(ntp2fit(ng+1),...
                                            dataList(group2List(ng+1,2),2)-dataList(group2List(ng+1,1),1)+1);
                                        
                                        %we do not check for minLength: we have
                                        %to do it at the very end
                                        
                                    else
                                        %there could be a new chunk. nextGroup
                                        %does not change
                                        
                                        dlEndIdx = group2List(ng+1,1)-1;
                                        
                                        isValid = (dlStartIdx < dlEndIdx);
                                        
                                    end
                                else
                                    %test valid below
                                end
                            end
                            
                            if isValid
                                %update lists:
                                %   group2List
                                %   ntp2fit
                                group2List(end+1,:) = [dlStartIdx,dlEndIdx];
                                ntp2fit(end+1,1)    = min([dataList(dlEndIdx,2)-dataList(dlStartIdx,1)+1],...
                                    ntp2fit(ng));
                            end
                        end %for chunk = 1:size(remainSE,1)
                    end %if isempty(remainSE)
                end %if isempty(winnerIdx)
                
            end % for ng = groups2fit' %check only for fitted groups!
            
            %remove entries
            %add all entries where ntp2fit <= minlenght
            rmGroupNumber = [rmGroupNumber;find(ntp2fit<=constants.MINLENGTH)];
            
            group2List(rmGroupNumber,:) = [];
            ntp2fit(rmGroupNumber)      = [];
            
            %update ngroup and restart loup
            nGroups = length(ntp2fit);
            
            if ~isempty(rmGroupNumber)
                %since we'll be looking at last/next groups: sort
                [group2List,sortIdx] = sortrows(group2List,[1,2]);
                ntp2fit = ntp2fit(sortIdx);
            end            
            %----END ACCEPT WINNERS
            
        end  %while nGroups > 0
        
        
        %END CASE 1
        
    case 2 %start with doublet fit
    otherwise %nothing here
end


%update dataListGroup(6:9)

%update only what we have fitted
idx2updateTime = find(dataListG(:,7)==0 );
idx2updateDist = find(dataListG(:,7)==0 & dataListG(:,10)==0); %no pause updates here!

if ~isempty(idx2updateTime)
    %update time
    dataListG(idx2updateTime,7) = time(dataListG(idx2updateTime,2),1)-time(dataListG(idx2updateTime,1),1);
end
if ~isempty(idx2updateDist)
    %update deltaDistance: calculate from fit: deltaD = slope*deltaT
    dataListG(idx2updateDist,9) = dataListG(idx2updateDist,4).*dataListG(idx2updateDist,7);
    %update delteDistanceSigma: dds = sqrt(deltaT)*slopeSigma
    dataListG(idx2updateDist,10) = sqrt(dataListG(idx2updateDist,7)).*max(dataListG(idx2updateDist,5:6),[],2);
end