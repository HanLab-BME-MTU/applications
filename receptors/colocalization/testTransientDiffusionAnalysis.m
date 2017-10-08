function [results,fract,switchCell] = testTransientDiffusionAnalysis(transDiff,modelDiff,timeStep)
%short script to evaluate accuracy of transient diffusion function
% Input:
% timestep normalizes length of model track to actual observed track
% Output:
% switchError-
% missError - 
% % 

%% Initialization
nanCount = 0;%Possibly remove
%Initialize outputs
fract = nan(length(transDiff),1);%fraction of correct track part
fract1Switch =  nan(length(transDiff),1);
switchCell = cell(length(transDiff),1); %Where switch happened, where it should have happened, if at all
testCoord = zeros(length(transDiff),5); %fraction of diffusion type present in track
modelCoord = zeros(length(transDiff),5);
noSwitch = 0;
ySwitch = 0;
ySwitchType = zeros(4,1);
noSwitchType = zeros(4,1);
switchDiff =NaN(length(transDiff),3); 
sd = 1; %switchDiff counter
errorBox = zeros(4,4);
prefBox = zeros(5,5,length(transDiff));
prefBoxNorm = zeros(5,5,length(transDiff));
prefBoxNorm1Switch = NaN(5,5,length(transDiff));
%% Analysis
for k = 1:length(transDiff)
    if isstruct(transDiff)
        track = transDiff(k).segmentClass.momentScalingSpectrum;
        lastT = size(track,1);%How many segments are there?
        % If not classified, then skip % REMOVING FOR NOW
    % %     if isnan(track(lastT,3))
    % %         nanCount = nanCount+1;
    % %         fract(k)=NaN;
    % %         continue
    % %     else
    % %     testT = zeros(track(lastT,2),1);
    % %     end

        %% Go through real data and replace track point with classification
        testT = zeros(track(lastT,2),1);    
        for j = 1:lastT
    %         if track(j,3) == 0 %Why is this removed?
    %             track(j,3) = 4;
    %         end
            testT(track(j,1):track(j,2),1) = track(j,3);
        end
    else
       testT = transDiff{k}';
       if ~isempty(testT(testT == 3))
           testT(testT==1) = 3;
           testT(testT==2) = 3;
%            testT(testT==4) = 0; %Work around to get correct results
       else
       testT(testT==2) = 4;
       testT(testT==1) = 2;
       testT(testT==4) = 0; %Work around to get correct results
       end
       testT(end+1) = testT(end); 
    end
    
    %% Go through ground truth and do the same
    model = modelDiff(k).segmentClass.momentScalingSpectrum;
    lastM = size(model,1);
    testM = zeros(model(lastM,2),1);
    for j = 1:lastM
        testM(model(j,1):model(j,2),1) = model(j,3);
    end
    testM = testM(1:(1/timeStep):end);
    %Skipping this part for now until use is found
    %
    %Change NaNs to some other value
    diffTypes = unique(testT(~isnan(testT)));% See how this works
    for d = 1: length(diffTypes)
        ind = find(testT(:,1)==diffTypes(d));
        testCoord(k,diffTypes(d)+1) = length(ind)/length(testT);
    end
    testN = isnan(testT(:,1));
    testT(testN) = 4;
    indT = find(testN);
    testCoord(k,5) = length(indT)/length(testT);
    
    
    diffTypes = unique(testM(~isnan(testM)));% See how this works
    for d = 1: length(diffTypes)
        ind = find(testM(:,1)==diffTypes(d));
        modelCoord(k,diffTypes(d)+1) = length(ind)/length(testM);
    end
    testN = isnan(testM(:,1));
    indM = find(testN);
    modelCoord(k,5) = length(indM)/length(testM);
    
    %}
%% General Error and error without no switch cases
    testFract = testM-testT;
    fract(k) = length(find(testFract==0))./length(testFract);
    %Old way: fract(k) = 1-sum(abs(modelCoord(k,:)-testCoord(k,:)))/2;
    if unique(testCoord(k,:))
    end
    fractFull = [testM,testT];
    typesD = unique(fractFull,'rows');
%     if find(typesD(:,2)==4)
%         typesD(find(typesD(:,2)==4),:) = [];
%     end
    for td = 1: size(typesD,1)
        prefBox(typesD(td,1)+1,typesD(td,2)+1,k)=length(find(ismember(fractFull,typesD(td,:),'rows')))./length(testM);
        prefBoxNorm(typesD(td,1)+1,typesD(td,2)+1,k)=length(find(ismember(fractFull,typesD(td,:),'rows')))./length(find(testM==typesD(td,1)));
    end
    %% Get types of error
    existSwitch = testT(2:end)-testT(1:end-1);
    sP = find(existSwitch~=0);
    switchTrack(k) = length(sP);
    if isempty(sP)
        switchPoint(:,1) = model(1,3); %Shows what was supposed to happen..
        if size(model,1)>1
        switchPoint(:,2)= model(2,3);
        else
         switchPoint(:,2)= model(1,3);   
        end
        switchPoint(1,3) = round(model(1,2).*timeStep);%---
        switchPoint(2,1) = testT(1);
        switchPoint(2,2) = testT(1);
        switchPoint(2,3) = NaN;
%         switchPoint(:,4) = track(1,3);
    else
        %% Get stats on single switch situations
        if length(sP) ==1
            fract1Switch(k) = fract(k);
            prefBoxNorm1Switch(:,:,k) = prefBoxNorm(:,:,k);
        end
        switchPoint(1,1) = model(1,3);
        if size(model,1)>1
        switchPoint(1,2) = model(2,3);
        else
        switchPoint(1,2) = model(1,3);  
        end
        switchPoint(1,3) = round(model(1,2).*timeStep);%---
        switchPoint(2:2+length(sP)-1,1) = testT(sP);% Monitor, changed from (sP-1), shouldn't make a difference
        switchPoint(2:2+length(sP)-1,2)= testT(sP+1);
        switchPoint(2:2+length(sP)-1,3) = sP;
%         switchPoint(:,4) = NaN;
    end
    if ~isnan(switchPoint(1,3)) && isnan(switchPoint(2,3)) %
        noSwitch = noSwitch+1;
        noSwitchType(switchPoint(2,2)+1,1) = noSwitchType(switchPoint(2,2)+1,1)+1;
        errorBox(switchPoint(1,1)+1,switchPoint(1,2)+1) = errorBox(switchPoint(1,1)+1,switchPoint(1,2)+1)+1;
    elseif ~isnan(switchPoint(1,3)) && ~isnan(switchPoint(2,3))
        amtSwitch =  switchPoint(1,3)-switchPoint(2:end,3);
        %Get either single switch or switch closest to real switch
        if length(amtSwitch) == 1
            switchDiff(k,1) = amtSwitch;
            if amtSwitch<0
                switchDiff(k,2) = switchPoint(1,1);
            elseif amtSwitch>0
                switchDiff(k,2) = switchPoint(1,2);
            else
                switchDiff(k,2) = switchPoint(1,1); %Technically shouldn't have "preference"
            end
            switchDiff(k,3) = 1;
        else
            switchDiff(k,1) = min(abs(amtSwitch));
            if min(abs(amtSwitch)) < 0
                switchDiff(k,2) = switchPoint(1,1);
            elseif min(abs(amtSwitch)) > 0
                switchDiff(k,2) = switchPoint(1,2);
            else
                switchDiff(k,2) = NaN;
            end
            switchDiff(k,3) = 0;
            sd=sd+1;
        end
    end
    switchCell{k} = switchPoint(:,:);
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
 results.NoSwitch = noSwitch./length(transDiff);
 results.switchTrack = switchTrack;
%  if find(errorBox)
%     [r,c] = find(errorBox);
%     for t = 1:2 %length(r)
%         results.TypeNoSwitch(t,1) =r(t)-1;
%         results.TypeNoSwitch(t,2) =c(t)-1;
%         results.TypeNoSwitch(t,3) = errorBox(r(t),c(t));
        results.TypeNoSwitch =noSwitchType./length(transDiff);
% %         results.TypeNoSwitch(1,1) =typesC(1);
% %         results.TypeNoSwitch(1,2) =typesC(2);
% %         results.TypeNoSwitch(1,3) = errorBox(typesC(1)+1,typesC(2)+1);
% %         results.TypeNoSwitch(2,1) =typesC(2);
% %         results.TypeNoSwitch(2,2) =typesC(1);
% %         results.TypeNoSwitch(2,3) = errorBox(typesC(2)+1,typesC(1)+1);
%     end
%  else
%      results.TypeNoSwitch =[0,0,0;0,0,0];
%  end
 results.avgFractionCorrect = [median(fract),prctile(fract,[25 75])];
 results.FullFractionCorrect = fract;%^previously showed mean and std
%  results.avgFractionCorrect1Switch = [nanmedian(fract1Switch),prctile(fract1Switch,[25 75])];%previously showed mean and std
 results.switchDiff = switchDiff;%[mean(switchDiff),std(switchDiff)];
%  results.switchDiff1Switch = switchDiff1Switch;
 results.FullFraction = prefBox;%mean(prefBox,3);
 results.FullFractionNorm = prefBoxNorm;%mean(prefBoxNorm,3);
 results.SwitchStats = switchCell;
%  results.FullFractionNorm1Switch = nanmean(prefBoxNorm1Switch,3);
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