 function [ output_args ] = gcaPlotHMMTransOfACertainType(toPlot,varargin)
% gcaPlotHMMTransOfACertainType 
% 
%%Input check
ip = inputParser;

ip.addParameter('OutputDirectory',pwd);

% Filters 
ip.addParameter('minStateTime', 90); % (in s) time in new state after transition 
ip.addParameter('minPercentTimePR',0.85); 
ip.addParameter('maxPercentLengthOutliers',0.1); % don't allow for more than 
% 10% of the tracks to be an outlier (outlier's currently IDed via median filtering of trajectory).
ip.addParameter('minGradBetweenStates',0.3); % minimum gradient between states to be considered:  
% currently calculated as the difference between percent pause/retract for
% state2 - state 1 
ip.addParameter('plotElongState',true); 
% plots
ip.addParameter('plotEndState',true); 

ip.parse(varargin{:});


%% 
nGroups = 1; 
% start by plotting the MDS Coords 
valuesAll = arrayfun(@(i) vertcat(toPlot.MDSValues{i}{:}),1:10,'uniformoutput',0); 
valuesAll = vertcat(valuesAll{:}); 
figure
scatter(valuesAll(:,1),valuesAll(:,2),20,[0.5,0.5,0.5],'filled'); 
hold on 
cmap = jet(20);  
count = 1; 
countOutliers =0; 
for iGroup = 1:nGroups 

    %     if ~isempty(projIDs{iGroup})
    projList = toPlot.info.projList{iGroup};
    nProjs = size(projList,1);
    for iProj = 1:nProjs
        
        % load Elong State
        % load StateVect
        %             cProj = projIDs{iGroup}(iProj);
        
        load([projList{iProj,1} filesep 'GrowthConeAnalyzer' filesep 'movieData.mat']);
        
        load([MD.outputDirectory_ filesep 'SegmentationPackage/StepsToReconstruct/IV_veilStem_length/Channel_1'...
            filesep 'veilStem.mat']);
        
        % calculate the outliers
        percentOutliers = sum(vertcat(veilStem.flagOutlier))/length(veilStem(:));
        out.percentOutliers(iProj) = percentOutliers; 
        if percentOutliers < ip.Results.maxPercentLengthOutliers
            
            load([MD.outputDirectory_ filesep ...
                'SegmentationPackage/StepsToReconstructTestBugFix20160426/GCAMeasurementExtraction_test20160510/WholeNeurite/' ...
                'Partition_Outgrowth_Trajectory_WithGrowthTimes_Spline0pt01'  filesep 'globalMeas.mat' ]);
            elongStates = globalMeas.outgrowth.stateAllFrames;
            elongStates = elongStates(1:119);
            name = gcaGetNeuriteID(MD.outputDirectory_);
            name = strrep(name,' ','_');
            
            load([MD.outputDirectory_ filesep 'HMMResults/20161025_InitialAnalysisNoVel' filesep ...
                'HMMResults_' name '.mat']);
            states = Results.ML_states';
            r = Results.ML_params.mu_emit;
            sigma = Results.ML_params.sigma_emit;
            % try first by replacing with percent pause
            nStates = length(unique(states));
            elongEachState =  arrayfun(@(x) elongStates(states == x),1:nStates,'uniformoutput',0);
            percentPauseRetract = cellfun(@(x) sum(x == 1 | x ==2)./length(x),elongEachState);
            %  percentPauseRetract = cellfun(@(x) sum(x ==2)./length(x),elongEachState);
            % find transitions where the difference between states is
            % larger than .7
            
            stateVectByPR = zeros(length(states),1);
            for i = 1:nStates
                stateVectByPR(states==i) = percentPauseRetract(i);
            end
            % find the indexing of all the state transitions
            transAll = diff(states);
            % state transitions are only those with a non-zero value, add
            % back the end point so can get time diff, (only care right now
            % about the time of the 2nd state so don't need 1). 
            % This will give the index of the last value before transition.
            % 
            transAllIdx = [find(transAll~=0);length(stateVectByPR)];
           
            % just get the PR calcs at the transition points. 
            stateVectByPRTrans= stateVectByPR(transAllIdx); 
            diffValues = diff(stateVectByPRTrans);
           % diffValuesAll = [find(diffValues~=0);length(stateVectByPR)];
%             idx = find(diffValues>ip.Results.minGradBetweenStates);
%             idx = [idx;length(stateVectByPR)];
           % get the state IDs for each transition in a cell 
            stateTransAll = arrayfun(@(x) states([transAllIdx(x) transAllIdx(x)+1]),...
                1:length(transAllIdx)-1,'uniformoutput',0);
          
                 %f length(idx)~=1
             % for now only care about the time of the state to which we
             % transition: so don't worry about counting the time of the first state 
             % Should have one time per transition.
             timeTrans = diff(transAllIdx);
             %timeTrans = arrayfun(@(x) diff(transAllIdx),1:length(transAllIdx)-1);
             % convert to sec for filtering
             timeTrans = timeTrans.*MD.timeInterval_;
             
              
                % screen the transitions. 
                for i = 1:numel(stateTransAll)
                    state1 = stateTransAll{i}(1);
                    state2 = stateTransAll{i}(2);
                    if (timeTrans(i) > ip.Results.minStateTime && percentPauseRetract(state2)>ip.Results.minPercentTimePR...
                            && diffValues(i)>ip.Results.minGradBetweenStates)
                        
                        scatter(r(1,state1),r(2,state1),50,cmap(iProj,:),'filled');
                        text(r(1,state1),r(2,state1),num2str(iProj),'color','k');
                        line([r(1,state1),r(1,state2)],[r(2,state1),r(2,state2)],'linewidth',1,'color',cmap(iProj,:));
                          
          
                        if ip.Results.plotEndState
                            gcaCircles(r(1,state2),r(2,state2),sigma(state2),'edgecolor','k','facecolor','none');
                            gcaCircles(r(1,state2),r(2,state2),2*sigma(state2),'edgecolor','k','facecolor','none');
                            
                        end
                        
                        %                     if ip.Results.plotEndStateScatterElong
                        %
                        %
                        %                     end
                        if ip.Results.plotElongState
                            MDSValuesC = toPlot.MDSValues{iGroup}{iProj};
                            valuesState = MDSValuesC(states==state2,:);
                            c(1) = 'c';
                            c(2) = 'b';
                            c(3) = 'r';
                            c(4) = 'm';
                            elongEndState = elongEachState{state2};
                            arrayfun(@(x) scatter(valuesState(elongEndState==x,1),valuesState(elongEndState==x,2),20,c(x),'filled'),1:4);
                        end
                        out.timeTrans(count) = timeTrans(i);
                        out.percentPR(count) = percentPauseRetract(state2);
                        out.projDir{count} = MD.outputDirectory_;
                        out.stateTrans{count} = stateTransAll{i}';
                        out.diffValues(count) = diffValues(i);
                        
                        
%                         e = elongStates(states==state2);  
%                         vector = [sum(e ==1)/length(e) sum(e==2)/length(e) sum(e==3)/length(e) sum(e==4)/length(e)] ; 
%                         out.percentState(count) = vector; 
                       
                        count = count+1;
                    %end
                    
                end
            end
            axis([-0.009,0.0025,-.0055,0.004]);
            axis('square');
              else
            countOutliers =1+countOutliers; 
            out.outliers{countOutliers} = MD.outputDirectory_; 
            
        end % if
      
        
    end % for iProjs
    %     end
end % for iGroups

saveas(gcf,[ip.Results.OutputDirectory filesep 'FilteredTransitions.eps'],'psc2'); 
saveas(gcf,[ip.Results.OutputDirectory filesep 'FilteredTransitions.fig']); 

params = ip.Results; 
save([ip.Results.OutputDirectory filesep 'params.mat'],'params'); 
save([ip.Results.OutputDirectory filesep 'out.mat'],'out'); 
