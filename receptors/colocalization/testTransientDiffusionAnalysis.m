function [fract,modelCoord,testCoord,switchCell] = testTransientDiffusionAnalysis(transDiff,modelDiff,timeStep)
%short script to evaluate accuracy of transient diffusion function
% timestep normalizes length of model track to actual observed track
% Output:
% switchError-
% missError - 
% % model = zeros(101,1);

%% Real Data Test
nanCount = 0;
fract = nan(length(transDiff),1);
switchCell = cell(length(transDiff),1);
testCoord = nan(length(transDiff),3);
modelCoord = nan(length(transDiff),3);
% modelCoord = nan(length(transDiff),3);
% figure; hold on;
for k = 1:length(transDiff)
    track = transDiff(k).segmentClass.momentScalingSpectrum;
    lastT = size(track,1);
    % If not classified, then skip
    if isnan(track(lastT,3))
        nanCount = nanCount+1;
        fract(k)=NaN;
        continue
    else
    testT = zeros(track(lastT,2),1);
    end
    %Go through real data and replace track point with classification
    for j = 1:lastT
%         if track(j,3) == 0
%             track(j,3) = 4;
%         end
        testT(track(j,1):track(j,2),1) = track(j,3);
    end
    %Go through ground truth and do the same
    model = modelDiff(k).segmentClass.momentScalingSpectrum;
    lastM = size(model,1);
    testM = zeros(model(lastM,2),1);
    for j = 1:lastM
        testM(model(j,1):model(j,2),1) = model(j,3);
    end   
    idT1 = find(testT(:,1)==1);
    idT2 = find(testT(:,1)==2);
    idT3 = find(testT(:,1)==0);
    idM1 = find(testM(:,1)==1);
    idM2 = find(testM(:,1)==2);
    idM3 = find(testM(:,1)==0);
    testCoord(k,1) = length(idT1)/length(testT);%Adding division to these
    testCoord(k,2) = length(idT2)/length(testT);
    testCoord(k,3) = length(idT3)/length(testT);
    
    modelCoord(k,1) = length(idM1)/length(testM);
    modelCoord(k,2) = length(idM2)/length(testM); 
    modelCoord(k,3) = length(idM3)/length(testM);
    
%     c = testM-testT;
%     idx = find(c==0);
%     fract(k) = length(idx)/length(c);
%% General Error
    fract(k) = 1-sum(abs(modelCoord(k,:)-testCoord(k,:)))/2;
    %% Get types of error
    existSwitch = testT(2:end)-testT(1:end-1);
    sP = find(existSwitch~=0);
    if isempty(sP)
        switchPoint(:,1) = model(1,3); %Shows what was supposed to happen..
        if size(model,1)>1
        switchPoint(:,2)= model(2,3);
        else
         switchPoint(:,2)= model(1,3);   
        end
        switchPoint(1,3) = model(1,2);
        switchPoint(2,1) = testT(1);
        switchPoint(2,2) = testT(1);
        switchPoint(2,3) = NaN;
%         switchPoint(:,4) = track(1,3);
    else
        switchPoint(1,1) = model(1,3);
        if size(model,1)>1
        switchPoint(1,2) = model(2,3);
        else
        switchPoint(1,2) = model(1,3)*timeStep;  
        end
        switchPoint(1,3) = model(1,2);
        switchPoint(2:2+length(sP)-1,1) = testT(sP-1);
        switchPoint(2:2+length(sP)-1,2)= testT(sP+1);
        switchPoint(2:2+length(sP)-1,3) = sP;
%         switchPoint(:,4) = NaN;
    end
    switchCell{k} = switchPoint;
    clear switchPoint
%      for j =1:track(lastT,2)
%          if c(j)==0
%              if testM(j) ==2
%                 plot(j-1:j,[k,k],'Color',[0 1 1]);
%              elseif testM(j) ==1
%                  plot(j-1:j,[k,k],'Color',[0 0 1]);
%              else
%                  plot(j-1:j,[k,k],'Color',[0 0 0]);
%              end
%          elseif c(j) ==1
%            plot(j-1:j,[k,k],'Color',[1 0 0]);  
%          elseif c(j) ==-1
%            plot(j-1:j,[k,k],'Color',[0 1 0]);  
%          end
%      end
end
% figure; 
% optimalHistogram(fract);



%% synthetic test
% removed for now.
% data = zeros(101,1);
% switchError = nan(length(transDiff),2);
% figure;
% p=zeros(1,2);
% for run =  1:length(transDiff)
%     
% % %     diffPop = size(trajTransDiffClass(run).segmentClass.momentScalingSpectrum,1);
% % % 
% % %     for r = 1: diffPop
% % %         model(trajTransDiffClass(run).segmentClass.momentScalingSpectrum(r,1):trajTransDiffClass(run).segmentClass.momentScalingSpectrum(r,2),1)...
% % %             =trajTransDiffClass(run).segmentClass.momentScalingSpectrum(r,3);  
% % %     end
% 
%     diffPopRaw = size(transDiff(run).segmentClass.momentScalingSpectrum,1);
%     if diffPopRaw >1
%         x= modelDiff(run).segmentClass.momentScalingSpectrum(1,2)-transDiff(run).segmentClass.momentScalingSpectrum(1,2);
%         
% %         dcRatio = transDiff(run).segmentClass.momentScalingSpectrum(1,4)/transDiff(run).segmentClass.momentScalingSpectrum(2,4);
%         if transDiff(run).segmentClass.momentScalingSpectrum(1,3) == 1
%             scatter(x,run, 'MarkerEdgeColor',[1 0 0],'MarkerEdgeColor',[1 0 0]);
%             switchError(run,1)=x;
%         elseif transDiff(run).segmentClass.momentScalingSpectrum(1,3) == 2 
%            scatter(x,run, 'MarkerEdgeColor',[0 0 1],'MarkerEdgeColor',[0 0 1]);
%            switchError(run,2)=x;
%         else
%             scatter(x,run, 'MarkerEdgeColor',[0 0 0],'MarkerEdgeColor',[0 0 0]);
%         end
% % %         for r = 1: diffPopRaw
% % %             data(transDiff(run).segmentClass.momentScalingSpectrum(r,1):transDiff(run).segmentClass.momentScalingSpectrum(r,2),1)...
% % %                 =transDiff(run).segmentClass.momentScalingSpectrum(r,3);  
% % %         end
%     else
%         if transDiff(run).segmentClass.momentScalingSpectrum(1,3) == 1
%             scatter(rand(1),0, 'MarkerEdgeColor',[1 0 0],'MarkerEdgeColor',[1 0 0]);
%             p(1,1) = p(1,1)+1;
%         elseif transDiff(run).segmentClass.momentScalingSpectrum(1,3) == 2 
%            scatter(rand(1),0, 'MarkerEdgeColor',[0 0 1],'MarkerEdgeColor',[0 0 1]);
%            p(1,2) = p(1,2)+1;
%         else
%             scatter(rand(1),0, 'MarkerEdgeColor',[0 0 0],'MarkerEdgeColor',[0 0 0]);
%         end
%        
%        
%     end
%     hold on
%     
% end
% missError = p;

end