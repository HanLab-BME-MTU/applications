function out = calcMTDynamics_main(inputData,verbose,TESTPROBLOW,TESTPROBHIGH,analysisType)
%calculate MT dynamics
%
%input(future): inputData(opt): structure with a field anaDat for every anaDat
%               to be analyzed. allowing job-files as input will follow
%               saveDir(opt): directory where the data should be saved 
%                             (0 if it should not be saved)
%               verbose(opt): whether to display anything (1,0)
%
%mtdDat(1:end-1) : data for each movie analyzed
%
%fields:  .data  .distance [distance, sigma] in mum
%                .time [time, sigma] in seconds
%
%               stateList contains all the data for analysis (point might not 
%                         necessarily be needed outside the program)
%                         NOTE: all slopes are given in microns/minute!
%
%         .stateList   .point [#of point, state (1/-1/0), slope, slopeSigma]
%                      .single [#of group (binary), state, startData#, endData#, deltaT, deltaD, slope, slopeSigma, significance]
%             (currently not used) .group = [counter, state(1/2/-1/-2/0), startIdx#, endIdx#, deltaT, deltaD, speed1, speedSigma1, speed2, speedSigma2...
%                                           avg. time in undetermined state, sum of single counters, significance1, significance2]
%
%               plotData contains all the data necessary to generate the plots
%               (currently not saved)
%
%         .plotData  fields .point, .single, .group1, .group2 with subfields
%                       xData1,xData2,xData3, yData1,yData2,yData3
%
%         .info .dataProperties
%
%         .statistics   cell array, structured {name}{mean value}{std}
%                       distance;
%                       catFreq, resFreq, growthSpeed, shrinkageSpeed,
%                       transition rates, avg distance travelled %time
%                       undetermined
%                       for group only
%
%mtdDat(end) .statistics: same statistics, but for all the analyzed movies
%
%
% 06-03-jonas

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%-----------------start calculations
numAnaDat = length(inputData);

for ct = 1:numAnaDat
    %read the correct anaDat
    anaDat = inputData(ct).anaDat;
    
    %calculate data
    
    %calc data. time = [time, timeSigma]
    data.time = [cat(1,anaDat.time),cat(1,anaDat.sigmaTime)];
    
    %only G1-data for now (later: look at labelcolor)
    selectedTags = [1,2];
    
    %find distances
    distanceMatrix=cat(3,anaDat.distanceMatrix);
    distance = squeeze(distanceMatrix(selectedTags(1),selectedTags(2),:));
    
    %extract data which is to be passed down for the sigma calculation
    distanceVectorMatrixN = cat(4, anaDat.distanceVectorMatrixN);
    distanceVector = squeeze(distanceVectorMatrixN(selectedTags(1),selectedTags(2),:,:))' .* (distance*ones(1,3));
    %calculate sigma
    [distanceSigma] = adgui_calcPlotData_distanceSigma(anaDat,selectedTags,distanceVector,distance,0);
    
    data.distance = [distance, distanceSigma];
    
    %clear up
    clear distanceMatrix distance distanceVectorMatrixN distanceVector distanceSigma
    
    %-----------------------------------------------------------------------------
    %for every triplet in the data: do a least squares fit using the known variances
    %of the parameters (only distance for now). Then go and check, whether the
    %slopes are definitely positive, negative or undecided, and group the different
    %MT growth states
    
    %use lscov for now. [x,sigmaX] = lscov(A,B,Sigma)
    
    nTimePoints = length(anaDat);
    
    %plot data
    if verbose
        figure('Name',anaDat(1).info.name,'NumberTitle','off');
        plot(data.time(:,1),data.distance(:,1),'b-d');
        myErrorbar(data.time(:,1),data.distance(:,1),[data.time(:,2),data.distance(:,2)]);
        hold on;
    end
    %set test probability. Confidence level is 1-2(1-X) (two-sided test)
    testValue = tinv(TESTPROBLOW,1);
    
    %initialize storing variables
    stateList = struct('point',zeros(nTimePoints,4),'single',[],'group',[]);
    %     xydat = struct('xData1',[],'xData2',[],'xData3',[],'yData1',[],'yData2',[],'yData3',[]);
    %     plotData = struct('point',xydat,'single',xydat,'group1',xydat,'group2',xydat);
    directionList = ones(nTimePoints,3)*NaN;
    
    %entry for first and last point in stateList
    stateList.point([1,end],1) = [1;nTimePoints];
    
    %loop through triplets, leave out first and last
    for i = 2:nTimePoints-1;
        %prepare matrices
        A = [data.time(i-1:i+1,1),ones(3,1)];
        B = data.distance(i-1:i+1,1);
        S = diag(data.distance(i-1:i+1,2).^2);
        
        %do least squares
        [X, sigmaX] = lscov(A,B,S);
        sigmaX = sigmaX + 1e-10;
        
        %find significant slopes
        %H0: slope = 0, H1: slope ~= 0            
        if abs(X(1))/sigmaX(1) > testValue
            %significant slope
            if X(1)>0
                direction = 1;
                col = 'g--';
                
                %stateList: #of point, state, slope in mum/min, slopeSigma
                stateList.point(i,:) = [i,1,X(1)*60,sigmaX(1)*60];
                
                %                 %save plotData according to state; end with NaN to allow
                %                 %discontinuos plot
                %                 plotData.point.xData1 = [plotData.point.xData1; A(:,1); NaN];
                %                 plotData.point.yData1 = [plotData.point.yData1; A*X; NaN];
                
            else 
                direction = -1;
                col = 'r--';
                
                %stateList: #of point, state, slope in mum/min, slopeSigma
                stateList.point(i,:) = [i,-1,X(1)*60,sigmaX(1)*60];
                
                %                 %save plotData according to state; end with NaN to allow
                %                 %discontinuos plot
                %                 plotData.point.xData2 = [plotData.point.xData2; A(:,1); NaN];
                %                 plotData.point.yData2 = [plotData.point.yData2; A*X; NaN];
                
            end % if X(1)>0
            
        else
            direction = 0;
            col = 'k--';
            
            %stateList: #of point, state, slope in mum/min, slopeSigma
            stateList.point(i,:) = [i,0,X(1)*60,sigmaX(1)*60];
            
            %             %save plotData according to state; end with NaN to allow
            %             %discontinuos plot
            %             plotData.point.xData3 = [plotData.point.xData3; A(:,1); NaN];
            %             plotData.point.yData3 = [plotData.point.yData3; A*X; NaN];
            
            
        end % if abs(X(1))/sigmaX(1) > TESTPRO
        
        if verbose
            %plot resulting line
            plot(A(:,1),A*X,col);
        end
        
        
    end %for i = 2:nTimePoints-1;
    
    %find connected regions of growth/shrinkage
    deltaPrev = abs([diff(stateList.point(:,2),1);NaN]);
    deltaNext = abs([NaN;NaN;stateList.point(2:end-1,2)-stateList.point(3:end,2)]); %there's no diff(v,-1)
    
    %wherever abs(deltaPrev)>abs(slopes(:,3)), either this one's a startIdx (if
    %delta=1) or the next one (if delta = 2). Same applies for deltaNext
    
    startIdx = find( deltaPrev > abs(stateList.point(:,2)) & deltaPrev == 1 ); %direction [0,1...]
    startIdx = [startIdx; find( deltaPrev > abs(stateList.point(:,2)) & deltaPrev == 2 ) + 1]; %direction [(-1),1...]
    startIdx = sort(startIdx); % there are no overlapping regions!
    
    stopIdx = find( deltaNext > abs(stateList.point(:,2)) & deltaNext == 1 ); %direction [...,1,0]
    stopIdx = [stopIdx; find( deltaNext > abs(stateList.point(:,2)) & deltaNext == 2 ) - 1]; %direction [...-1,(1...)]
    stopIdx = sort(stopIdx); % there are no overlapping regions!
    
    %take care of border effects: if #2 is significant, it still does not appear in
    %start yet - include #1, too!
    if startIdx(1)>=stopIdx(1)
        startIdx = [1;startIdx];
    end
    %same for end
    if stopIdx(end)<=startIdx(end)
        stopIdx = [stopIdx;nTimePoints];
    end
    
    %make sure the startIdx/stopIdx contain all points
    if startIdx(1)~= 1
        stopIdx = [startIdx(1);stopIdx];
        startIdx = [2;startIdx];
    end
    if stopIdx(end) ~= nTimePoints
        startIdx = [startIdx;stopIdx(end)];
        stopIdx = [stopIdx;nTimePoints - 1];
    end
    
    
    %-------------------calculate statistics
    %loop through all start/stop pairs (start at #1 and end with #end) and the
    %gaps in between and classify
    
    idxList = [startIdx,stopIdx];
    idxListFull = []; %list containing all regions; will be passed to the group fitting
    %init counters
    growthCt = 1;
    shrinkageCt = 1;
    undeterminedCt = 1;
    
    
    
    %loop throug indexList, remove entries once they have been evaluated
    while ~isempty(idxList)
        
        %select current data region...
        currentBlock = [idxList(1,:)];
        
        
        %...and delete it from the list
        idxList(1,:) = [];
        
        %insert undetermined region if necessary
        if ~isempty(idxList);
            if currentBlock(2) < idxList(1,1)
                idxList = [currentBlock(2),idxList(1,1);idxList];
            end
        end
        
        %make sure this is no block with zero length
        if currentBlock(1) < currentBlock(2)
            
            %test if we want to merge with subsequent block: if end of
            %currentBlock is equal to the start of the next block and the
            %slopes point into the same direction: merge
            if ~isempty(idxList)
                while currentBlock(2)==idxList(1,1) & stateList.point(idxList(1,1)+1,2) == stateList.point(currentBlock(1)+1,2)
                    
                    %make block larger
                    currentBlock(2) = idxList(1,2);
                    
                    %delete entry from list
                    idxList(1,:) = [];
                    
                    %be sure to leave loop if idxList is now empty
                    if isempty(idxList)
                        break
                    else %insert undetermined region if necessary
                        if currentBlock(2) < idxList(1,1)
                            idxList = [currentBlock(2),idxList(1,1);idxList];
                        end
                    end
                    
                end
            end
            
            %add the current block to the new idxList for groupFitting
            idxListFull = [idxListFull;currentBlock];
            
            %lengthCurrentBlock: # of data points, not # of segments!
            lengthCurrentBlock = diff(currentBlock)+1;
            
            %classify according to direction of startIdx+1 (at startIdx, the state is generally undetermined).
            currentState = stateList.point(currentBlock(1)+1,2);
            switch abs(currentState)*(lengthCurrentBlock>2) %we need at least 3 data points for linear fit!
                
                case 1 %growth/shrinkage
                    
                    %select growth or shrinkage and get specific data
                    switch currentState
                        case 1 %growth
                            col = '-g';
                            counterStr = 'growthCt';
                        case -1 %shrinkage
                            col = '-r';
                            counterStr = 'shrinkageCt';
                    end
                    
                    
                    %calculate slope:
                    %prepare matrices
                    A = [data.time(currentBlock(1):currentBlock(2),1),ones(lengthCurrentBlock,1)];
                    B = data.distance(currentBlock(1):currentBlock(2),1);
                    S = diag(data.distance(currentBlock(1):currentBlock(2),2).^2);
                    
                    %do weighted least squares
                    [X, sigmaX] = lscov(A,B,S);
                    
                    %add something small to avoid division by zero
                    sigmaX = sigmaX + 1e-10;
                    
                    %test sigmificance. 
                    degF = lengthCurrentBlock - 2;
                    testValue = tinv(TESTPROBHIGH,degF);
                    
                    if abs(X(1))/sigmaX(1)>testValue & ...
                            imag(sigmaX(1))==0 & real(sigmaX(1)~=1e-10)
                        sigmaX = sigmaX - 1e-10;
                        if verbose
                            %plot result
                            plot(A(:,1),A*X,col);
                        end
                        %count only if significant
                        eval(['counter = ',counterStr,';']);
                        
                        %update counter
                        eval([counterStr,'=',counterStr,'*2;']);
                        
                    else
                        sigmaX = sigmaX - 1e-10;
                        if verbose
                            %plot result
                            plot(A(:,1),A*X,':k');
                        end
                        %count as not significant
                        counter = undeterminedCt;
                        undeterminedCt = undeterminedCt * 2;
                        currentState = 0;
                        
                    end
                    
                    %transform speeds into microns/minutes
                    speedMuM = 60*X(1);
                    speedSigmaMuM = 60*(sigmaX(1));
                    
                    %                     %store results
                    %                     plotData.single.xData1 = [plotData.single.xData1; A(:,1); NaN];
                    %                     plotData.single.yData1 = [plotData.single.yData1; A*X; NaN];
                    
                    %stateList: #of event, state, startIdx#, endIdx#, deltaT,
                    %deltaD, speed, speedSigma, significanceSwitch
                    deltaT = A(end,1) - A(1,1);
                    deltaD = B(end,1) - B(1,1);
                    stateList.single(end+1,:) = [counter, currentState, currentBlock, deltaT, deltaD, speedMuM, speedSigmaMuM];
                    
                case 0 %undetermined
                    
                    %no calc slope
                    A = [data.time(currentBlock(1):currentBlock(2),1)];
                    B = data.distance(currentBlock(1):currentBlock(2),1);
                    
                    %store results
                    %plotData.single.xData3 = [];
                    %plotData.single.yData3 = [];
                    
                    %stateList: #of event, state, startIdx#, endIdx#, deltaT,
                    %deltaD, speed, speedSigma, significanceSwitch
                    deltaT = A(end,1) - A(1,1);
                    deltaD = B(end,1) - B(1,1);
                    stateList.single(end+1,:) = [undeterminedCt, 0, currentBlock, deltaT, deltaD, 0, 0];
                    
                    %update counter
                    undeterminedCt = undeterminedCt * 2;
                    
            end %switch stateList.point(currentBlock(1)+1,2)
            
        end %if currentBlock(1) < currentBlock(2)
        
    end %while ~isempty(idxlist) #1
    
    
    
    
%     %===========================================================================
%     
%     %now loop through the full idxlist: if there are regions of e.g.
%     %growth-undet-growth, consider them as one partially interrupted growth.
%     %Do two linefits: one with equal slope but different intercept for all, and
%     %one with both equal slope and intercept
%     
%     
%     %init counters
%     growthCt = 1;
%     shrinkageCt = 1;
%     undeterminedCt = 1;
%     
%     %loop through it
%     while ~isempty(idxListFull)
%         
%         %select current data region...
%         currentBlock = [idxListFull(1,:)];
%         
%         %...and delete it from the list
%         idxListFull(1,:) = [];
%         
%         %no need to insert undetermined region: all are in there already
%         
%         %make sure this is no block with zero length
%         if currentBlock(1) < currentBlock(2)
%             
%             %get current state
%             currentState = stateList.point(currentBlock(1)+1,2);
%             
%             %no need to test merging: this has been done already
%             
%             %now test if next block is the same direction. No testing yet
%             %whether they are separated by more than a minimum distance
%             
%             %reset counters
%             currentGroupIdx = find(stateList.single(:,3)==currentBlock(1));
%             if ~isempty(currentGroupIdx)
%                 totalCt = stateList.single(currentGroupIdx,1);
%                 undefCt = 0.5;
%             else
%                 totalCt = 0.5;
%                 undefCt = 0.25;
%             end
%             undetLength = 0;
%             
%             if ~isempty(idxListFull) & currentState ~= 0
%                 %add blocks while either the next block or the next but one have
%                 %the same state (and the one in between is undetermined and not longer than two intervals
%                 nextGood = (stateList.point(idxListFull(1,1)+1,2) == currentState);
%                 if size(idxListFull,1)>1
%                     nextButOneGood = stateList.point(idxListFull(1,1)+1,2) == 0 & stateList.point(idxListFull(2,1)+1,2) == currentState ...
%                         & idxListFull(2,1) - currentBlock(2) < 3;
%                 else
%                     nextButOneGood = 0;
%                 end
%                 
%                 while nextGood | nextButOneGood
%                     
%                     %add block
%                     currentBlock = [currentBlock;idxListFull(1,:)];
%                     
%                     %delete entry from list
%                     idxListFull(1,:) = [];
%                     
%                     %remember contributing block (add number if exist (i.e. was long enough to be fitted), else add fraction)
%                     %but only if it is a determined state
%                     if stateList.point(currentBlock(end,2)-1,2) ~= 0
%                         currentGroupIdx = find(stateList.single(:,3) == currentBlock(end,1));
%                         if ~isempty(currentGroupIdx)
%                             totalCt = totalCt + stateList.single(currentGroupIdx,1);
%                         else
%                             totalCt = totalCt + undefCt;
%                             undefCt = undefCt/2;
%                         end
%                     else
%                         %remember only the lenght of the undetermined states
%                         undetLength = undetLength + diff(currentBlock(end,:));
%                     end
%                     
%                     %be sure to leave loop if idxList is now empty
%                     if isempty(idxListFull)
%                         break
%                     end %if stateList.point(currentBlock(end,2)-1,2) ~= 0
%                     
%                     %check next blocks
%                     nextGood = (stateList.point(idxListFull(1,1)+1,2) == currentState);
%                     if size(idxListFull,1)>1
%                         nextButOneGood = stateList.point(idxListFull(1,1)+1,2) == 0 & stateList.point(idxListFull(2,1)+1,2) == currentState ...
%                             & idxListFull(2,1) - currentBlock(2) < 3;
%                     else
%                         nextButOneGood = 0;
%                     end
%                     
%                 end %while nextGood | nextButOneGood
%                 
%             end %if ~isempty(idxListFull) & currentState ~= 0
%             
%             
%             %-----------------------------
%             %if we have found a group: do grouped line fit. Else, just use the
%             %results from above
%             switch abs(currentState)*(currentBlock(end)-currentBlock(1)>2)
%                 case 1 %growth or shrinkage
%                     
%                     %select growth or shrinkage and get specific data
%                     switch currentState
%                         case 1 %growth
%                             col = 'g';
%                             counterStr = 'growthCt';
%                         case -1 %shrinkage
%                             col = 'r';
%                             counterStr = 'shrinkageCt';
%                     end
%                     
%                     
%                     
%                     %check whether it is a group
%                     if size(currentBlock,1) == 1
%                         %use previous fit
%                         sLIdx = find(stateList.single(:,3)==currentBlock(1));
%                         
%                         if stateList.single(sLIdx,2) ~= 0
%                             %count only if significant
%                             eval(['counter = ',counterStr,';']);
%                             %update counter
%                             eval([counterStr,'=',counterStr,'*2;']);
%                         else
%                             counter = undeterminedCt;
%                             undeterminedCt = undeterminedCt * 2;
%                         end
%                         
%                         %stateList.group = [counter, state, startIdx#, endIdx#, deltaT, deltaD, speed1, speedSigma1, speed2, speedSigma2...
%                         %                       time in undetermined state, sum of single counters, significance1, significance2]
%                         stateList.group(end+1,:) = [counter,stateList.single(sLIdx,2:8),0,0,0,stateList.single(sLIdx,[1]),abs(stateList.single(sLIdx,2)),0];
%                         
%                         %store plotData
%                         %                         plotData.group1.xData1 = [plotData.group1.xData1; A1(:,1); NaN];
%                         %                         plotData.group1.yData1 = [plotData.group1.yData1; A1*X1; NaN];
%                         
%                         
%                         
%                     else
%                         %calculate grouped fit. Every even entry in currentBlock
%                         %is a undetermined phase, while every odd entry is a
%                         %growth phase
%                         
%                         nEntries = size(currentBlock,1);
%                         nOddEntries = ceil(nEntries/2);
%                         
%                         %calculate the time spent in undetermined state
%                         deltaTu = 0;
%                         for i = nEntries-1:-2:2
%                             deltaTu = deltaTu + data.time(currentBlock(i,2),1) - data.time(currentBlock(i,1),1);
%                             %remove undetermined entries
%                             currentBlock(i,:) = [];
%                         end
%                         %deltaTuAvg = deltaTu/(nOddEntries-1);
%                         
%                         %prepare fitting
%                         subBlockLength = diff(currentBlock,1,2)+1;
%                         numDataPoints = sum(subBlockLength);
%                         A1 = zeros(numDataPoints,nOddEntries+1);
%                         A2 = ones(numDataPoints,2);
%                         B = [];
%                         S = [];
%                         lastIdx = 0;
%                         
%                         for i = 1:nOddEntries
%                             %set time data and ones in the right column (we allow different intercepts here!)
%                             A1(lastIdx+1 : lastIdx+subBlockLength(i), [1,i+1]) = [data.time(currentBlock(i,1):currentBlock(i,2),1),ones(subBlockLength(i),1)];
%                             B = [B;data.distance(currentBlock(i,1):currentBlock(i,2),1)];
%                             S = blkdiag(S,diag(data.distance(currentBlock(i,1):currentBlock(i,2),2).^2));
%                             lastIdx = lastIdx + subBlockLength(i);
%                         end
%                         
%                         %write time data for A2
%                         A2(:,1) = A1(:,1);
%                         
%                         %fit lines
%                         [X1,sigmaX1] = lscov(A1,B,S);
%                         [X2,sigmaX2] = lscov(A2,B,S);
%                         sigmaX1 = sigmaX1 + 1e-10;
%                         sigmaX2 = sigmaX2 + 1e-10;
%                         
%                         %test sigmificance. 
%                         degF1 = numDataPoints - 1 - nOddEntries;
%                         degF2 = numDataPoints - 1 - 1; %only 1 slope fitted!
%                         testValue1 = tinv(TESTPROBHIGH,degF1);
%                         testValue2 = tinv(TESTPROBHIGH,degF2);
%                         
%                         if abs(X1(1))/sigmaX1(1)>testValue1
%                             significance1 = 1;
%                             if verbose
%                                 %plot result
%                                 %plot(A1(:,1),A1*X1,[col,'+-']);
%                             end
%                             %count only if significant
%                             eval(['counter = ',counterStr,';']);
%                             %update counter
%                             eval([counterStr,'=',counterStr,'*2;']);
%                             
%                             
%                         else
%                             significance1 = 0;
%                             if verbose
%                                 %plot result
%                                 %plot(A1(:,1),A1*X1,'k+-');
%                             end
%                             counter = undeterminedCt;
%                             undeterminedCt = undeterminedCt * 2;
%                             currentState = 0;
%                             
%                         end
%                         
%                         %X2 does not yet decide state
%                         if abs(X2(1))/sigmaX2(1)>testValue2
%                             significance2 = 1;
%                             if verbose
%                                 %plot result
%                                 %plot(A2(:,1),A2*X2,[col,'o-']);
%                             end
%                         else
%                             significance2 = 0;
%                             if verbose
%                                 %plot result
%                                 %plot(A2(:,1),A2*X2,'ko-');
%                             end
%                         end
%                         
%                         %transform speeds into microns/minutes
%                         speedMuM1 = 60*X1(1);
%                         speedSigmaMuM1 = 60*(sigmaX1(1));
%                         speedMuM2 = 60*X2(1);
%                         speedSigmaMuM2 = 60*(sigmaX2(1));
%                         
%                         %                         %store results
%                         %                         plotData.group1.xData1 = [plotData.group1.xData1; A1(:,1); NaN];
%                         %                         plotData.group1.yData1 = [plotData.group1.yData1; A1*X1; NaN];
%                         %                         plotData.group2.xData1 = [plotData.group2.xData1; A2(:,1); NaN];
%                         %                         plotData.group2.yData1 = [plotData.group2.yData1; A2*X2; NaN];
%                         
%                         deltaT = A1(end,1) - A1(1,1);
%                         deltaD = B(end,1) - B(1,1);
%                         
%                         %stateList.group = [counter, state, startIdx#, endIdx#, deltaT, deltaD, speed1, speedSigma1, speed2, speedSigma2...
%                         %                       time in undetermined state, sum of single counters, significance1, significance2]
%                         stateList.group(end+1,:) = [counter, 2*currentState, currentBlock(1), currentBlock(end), deltaT, deltaD,...
%                                 speedMuM1, speedSigmaMuM1, speedMuM2, speedSigmaMuM2, deltaTu, totalCt, significance1, significance2];
%                         
%                         
%                     end %size(currentBlock,1) == 1
%                     
%                     
%                     
%                 case 0 %undetermined
%                     %always use previous fit
%                     sLIdx = find(stateList.single(:,3)==currentBlock(1));
%                     
%                     %stateList.group = [counter, state, startIdx#, endIdx#, deltaT, deltaD, speed1, speedSigma1, speed2, speedSigma2...
%                     %                       time in undetermined state, sum of single counters, significance1, significance2]
%                     stateList.group(end+1,:) = [undeterminedCt,stateList.single(sLIdx,2:8),0,0,0,stateList.single(sLIdx,[1]),abs(stateList.single(sLIdx,2)),0];
%                     
%                     %store plotData
%                     %                         plotData.group1.xData3 = [plotData.group1.xData3; A1(:,1); NaN];
%                     %                         plotData.group1.yData3 = [plotData.group1.yData3; A1*X1; NaN];
%                     
%                     %update growthCt
%                     undeterminedCt = undeterminedCt *2;
%             end
%             
%         end %if currentBlock(1) < currentBlock(2)
%         
%     end %while ~isempty(idxListFull) #1
%     
%     
    
    %----------------------calculate overall movie statistics-------------------------
    %for now, use the stateList.single information (but don't discard the rest!)
    
    %test for values that are too good to be true

    %use log, because we want to see differences in the order of magnitude
    sigmaList(1).idx = find(stateList.single(:,2)>0);
    sigmaList(1).sigmaL = log(stateList.single(sigmaList(1).idx,8));
    sigmaList(2).idx = find(stateList.single(:,2)<0);
    sigmaList(2).sigmaL = log(stateList.single(sigmaList(2).idx,8));
%     sigmaList(3).idx = find(stateList.single(:,2)==0);
%     sigmaList(3).sigmaL = log(stateList.single(sigmaList(3).idx,8));
    
    
    %loop through the sigmaLists and do a least median squares outlier
    %detection
    for i = 1:length(sigmaList)
        if ~isempty(sigmaList(i).idx)
            %calculate median
            medAll = median(sigmaList(i).sigmaL);
            deltaSq = (sigmaList(i).sigmaL-medAll).^2;
            medDeltaSq = median(deltaSq);
            %according to Rousseeuw & Leroy, 1987, comparing the squared
            %differences from the median to the median times a magic number
            %shows which values can be considered outliers
            magicNumber2=1.4826^2;
            testValue = deltaSq/(magicNumber2*medDeltaSq);
            badIdx = find(testValue>=9);
            %write NaNs instead of -1/0/+1 into stateList.single
            stateList.single(sigmaList(i).idx(badIdx),2) = NaN;
        end
    end
    
    [statistics,statisticsCell] = calcMTDynamics_movieStats(data.distance,stateList.single,verbose,anaDat(1).info.name);
    
    
    %store date in main structure
    
    mtdDat(ct).statistics = statistics;
    mtdDat(ct).statisticsCell = statisticsCell;
 
    mtdDat(ct).stateListSingle = stateList.single; %do this separately (we want to use cat later)
    mtdDat(ct).stateList.point = stateList.point;
    mtdDat(ct).stateList.single = stateList.single;
%     mtdDat(ct).stateList.group = stateList.group;
    mtdDat(ct).data = data;
    mtdDat(ct).info.dataProperties = anaDat(1).info.dataProperties;
    mtdDat(ct).info.name = anaDat(1).info.name;
    mtdDat(ct).info.probabilities = [TESTPROBLOW, TESTPROBHIGH];
    
end %for ct = 1:numAnaDat


%-----------------calculate overall statistics----------------------------

%decide what to do according to analysisType
%init out
out = {};
numDatMinus = 0;

for type = bsum2bvec(analysisType)'
    
    switch type
        case 1 
            %calc mtdDat
            
            %collect data of individual trajectories: Use individual measurements for
            %global measurement (i.e. make one huge trajectory)
            
            %collect distance data
            numData = length(mtdDat);
            dataAll = [];
            for i = 1:numData
                dataAll = [dataAll;mtdDat(i).data.distance];
            end
            %collect stateList data
            stateListAll = cat(1,mtdDat.stateListSingle);
            
            [statisticsAll,statisticsAllCell] = calcMTDynamics_movieStats(dataAll,stateListAll,verbose,'all files');
            
            mtdDat(end+1).statistics = statisticsAll;
            mtdDat(end).statisticsCell = statisticsAllCell;
            out = [out,{mtdDat}];
            
            %there is one more data entry now
            numDatMinus = 1;
            
        case 2
            %calc conDat
            
            %loop through individual trajectory data: do first overall stats on the
            %first, then on the first and second trajectory etc
            
            %careful! mtdDat could be longer!
            numData = length(mtdDat)-numDatMinus;
            dataAll = [];
            
            for i = 1:numData
                %collect data
                dataAll = [dataAll;mtdDat(i).data.distance];
                stateListAll = cat(1,mtdDat(1:i).stateListSingle);
                
                %calculate statistics
                statisticsAll = calcMTDynamics_movieStats(dataAll,stateListAll,verbose,[num2str(i),' files']);
                
                %assign conDat
                conDat(i).statistics = statisticsAll;
                
            end
            
            out = [out,{conDat}];
            
    end %switch type
    
end %for type = bsum2bvec(analysisType)'


















