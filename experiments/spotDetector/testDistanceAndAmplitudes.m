function [nParms,nQAll,flag,rmIdx] = testDistanceAndAmplitudes(parms,QAll,chi,dataProperties,isNplus1,degreesOfFreedom)
% TESTDISTANCE tests whether the distances between the spots and their amplitudes are significant
%
% SYNOPSIS [nParms,nQAll,flag,rmIdx] = testDistanceAndAmplitudes(coord,QAll,chi,amp,dataProperties,isNplus1)
%
% INPUT   parms   : list of all fitparameters
%         QAll    : Q matrix of fitted model
%         chi     : chi squared of fitted model (=sigmaZeroHat)
%         dataProperties: constant definitions for project
%         isNplus1: whether this is a N+1-fit (the new spot would be in the first position)
%         degreesOfFreedom: degrees of freedom of the fit
%
% OUTPUT nParms : significant parameters
%        n*     : other parameters of significant coordinates
%        flag   : indices of deleted spots
%        rmIdx  : indices of entries in ub,lb to remove
%
%c: 11/03 jonas (based on testdistance by dT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%CONST DEFINITIONS
T_TEST_PROB=dataProperties.T_TEST_PROB;

%init parameters
nSpots=floor(length(parms)/4); %number of spots (last parm is bg)
endList = 4*nSpots;

flag = []; %init flag

%read amps and coords
coordIdx = sort([1:4:endList 2:4:endList 3:4:endList]);
ampIdx = [4:4:endList];

coord = reshape(parms(coordIdx),[length(coordIdx)/3,3]);
amp   = parms( ampIdx );

%chop Q into Qc (coords) and Qa (amps)
Qc = QAll(coordIdx,coordIdx);
Qa = QAll(ampIdx,ampIdx);

%---------FIRST TEST: test every spot in nCoord for zero amplitude
for i = 1:nSpots
    %null-hypothesis: the amplitude = 0
    %H1: amplitude ~= 0
    %testValue = amplitude/sqrt(Qa*sigmaZeroHat)
    testValue = amp(i)/sqrt(Qa(i,i)*chi);
    
    disp(sprintf('%f, %f, %f', testValue, tinv(1-(T_TEST_PROB),1),tinv(1-(T_TEST_PROB),76)))
    
    %if H0 accepted, we throw the spot out (we don't want spots with zero amplitude!)
    if testValue < tinv(1-(T_TEST_PROB),degreesOfFreedom) %one-sided test (amp is not going to be < 0)
        flag = [flag;i];
        snr = amp(i)/sqrt(chi);
    else
        %all is good
    end
end
%-----END FIRST TEST: test every spot in nCoord for zero amplitude


%---------2nd TEST: test every pair for zero distance
% crosscheck every pair for zero distance: if there are two spots that we
% can not distinguish from one another, remove the spot that has the bigger
% Q-matrix. If there are three spots and only one pair fails, remove the
% worse of the two only if there has no single amp been deleted.
% of course, if a spot has been added through a N+1-fit, it is deleted
% faster than others

%store the indistinguishable pairs
indistinguishable = [];
sigD = repmat(NaN,[nSpots,nSpots]);

for i=1:(nSpots-1) %test 1 against 2,3,4; 2 against 3,4; 3 against 4
    if ~any(i == flag) %don't test a spot that has been deleted already
        for j=(i+1):nSpots
            if ~any(j == flag) %don't test a spot that has been deleted already
                
                %read corresponding Q-matrix entries
                Qxx=blkdiag(Qc((i-1)*3+1:(i-1)*3+3,(i-1)*3+1:(i-1)*3+3), Qc((j-1)*3+1:(j-1)*3+3,(j-1)*3+1:(j-1)*3+3));
                
                %do the standard gaussian error propagation to get Qdd
                difference=(coord(i,:)-coord(j,:));
                distance=sqrt(sum(difference.^2));
                H=1/distance*[difference -difference];
                Qdd=H*Qxx*H';
                chiT=chi; %all jointly fitted spots have the same chi
                
                %calculate the denominator of the testValue
                sigD(i,j)=sqrt(Qdd*chiT);
                
                testValue=distance/sigD(i,j);
                
                %test: H0: zero distance
                %H1: nonzero distance (positive->one-sided test)
                %keep only coords that pass as nonzero distance
                
                % degrees of freedom stays constant even though we delete
                % individual spots: The df for the test are the ones used
                % in the fitting
                if testValue<tinv(1-(T_TEST_PROB),degreesOfFreedom) %one-sided
                    indistinguishable = [indistinguishable,[i;j]];
                end
            end
        end
    end
end

%now we have to decide which spots to remove
if ~isempty(indistinguishable)
    %count which spot appears the most
    [uniqueEntries,numberOfOccurences] = countEntries(indistinguishable);
    spotCount = [uniqueEntries,numberOfOccurences];
    spotCount = -sortrows(-spotCount,2); %sort in descending order
    
    %if the maximum number of entries is smaller than nsp-1 and we have
    %already removed a spot above, don't remove yet.
    %else, remove the one with the most entries. If there are two or more spots,
    %the N+1 of the two or the one with the largest Q-matrix is removed
    flagLength = length(flag);
    nEqual = length(find(spotCount(:,2)==spotCount(1,2)));
    switch 1*(flagLength~=0) + 2*(spotCount(1,2)==nSpots-flagLength) + 4*(nEqual>1)
        %case 0: no spots have been removed & less than max indist. & only one
        %max -> remove that max
        %case 1: spots have been removed & less than max indist. & only one
        %max -> do nothing
        %case 2: no spots have been removed & max indist. & only one max ->
        %remove
        %case 3: spots have been removed & max indist. & only one max ->
        %remove
        %case 4+ check & remove except 5
        
        case {1,5}
            %don't remove anything
            
        case {0,2,3}
            %remove the only max
            flag = [flag;spotCount(1,1)];
            
        case {4,6,7}
            %remove one of many: check which
            %1) any n+1?
            %2) any maxQ?
            
            %1) n plus 1 is the first entry in coordList
            if isNplus1 & any(spotCount(1:nEqual,1) == 1)
                flag = [flag;1];
            else
                %loop to read Qc traces
                for i = 1:nEqual
                    QTrace(i) = trace(Qc(3*(i-1)+1:3*i,3*(i-1)+1:3*i));
                end
                [dummy,maxIdx] = max(QTrace);
                flag = [flag;spotCount(maxIdx(1),1)]; %the third of the first three is number three. the (1) is just in case...
            end
        otherwise
            %there should be nothing that comes here
    end %switch
end %if ~isempty(indistinguishable)


%if we have a nonempty flag: we remove the spot to return nParms,nQAll
rmIdx = [];
if ~isempty(flag)
    %find entries in goodIdx we have to remove
    for i = 1:length(flag)
        rmIdx = [rmIdx,4*(flag(i)-1)+1:4*flag(i)];
    end
end

%assign output and remove what we have to
nParms = parms;
nQAll = QAll;
nParms(rmIdx) = [];
%rm in two steps, because with Q(rm,rm) we would cut a hole into the matrix
nQAll(rmIdx,:) = [];
nQAll(:,rmIdx) = [];
