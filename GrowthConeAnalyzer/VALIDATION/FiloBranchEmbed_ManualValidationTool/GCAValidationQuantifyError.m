function [ output_args ] = GCAValidationQuantifyError(reconstructDir,varargin)
%GCAValidationCalculatePercentError 
%% Runs through the manual output and calculates a number of error stats
% output is a heat map 

% note(used to be GCAValidationCalculatePercentErrorWithBranch)
%% 
ip = inputParser;

ip.CaseSensitive = false;
ip.addRequired('reconstructDir');
defaultType{1} = 'VeilFilo';
defaultType{2} = 'Branch';
defaultType{3} = 'Embed';
ip.addParameter('type',defaultType);
ip.addParameter('OutputDirectory',[]); 

ip.parse(reconstructDir,varargin{:});
%%
typesAll = ip.Results.type;

if isempty(ip.Results.OutputDirectory)
    outPutDirectory = pwd; 
else 
    outPutDirectory = ip.Results.OutputDirectory; 
end 

x = dir(reconstructDir);
x = x(3:end);
list = arrayfun(@(i) x(i).name,1:length(x),'uniformoutput',0);

% will eventually have a filter here that will check if done and only load
% those not finished if user desires.


idxInclude  = listSelectGUI(list,[],'move');


% listSelectGUI to select a file to work on
toTest = list(idxInclude);

groupIDs = cellfun(@(x) getGroupIDs(x),toTest,'uniformoutput',0)';  
% order{1} = 'KDNo';
% order{2} = 'KDCDC';
% order{3} = 'KDRhoA';
% order{4} = 'KDsrGAP';
% order{5} = 'KDBetaPix';
% order{6} = 'KDTrio';
% order{7} = 'KDDock';
% order{8} = 'CAAX';
sampled = unique(groupIDs,'stable'); 

count = 1; 
for iType = 1:numel(typesAll)
    typeC = typesAll{iType}; 
    collectMat = zeros(numel(toTest),5); 
for iProj = 1: numel(toTest)
    
    %if  strcmpi(ip.Results.type,'veilFilo');
    load([reconstructDir filesep toTest{iProj} filesep ...
        'Overlay' typeC filesep 'filoInfo.mat']);
    FN = load([reconstructDir filesep toTest{iProj} filesep ...
        'Overlay' typeC filesep 'ValidationOverlays' filesep ...
        'False_Neg_' typeC '_Coords.mat']);
    FP = load([reconstructDir filesep toTest{iProj} filesep ...
        'Overlay' typeC filesep 'ValidationOverlays' filesep ...
        'False_Pos_' typeC '_Coords.mat']);
    
    
    
    % just collect for now
    collectMat(iProj,1) = length(filoInfoFilter); % filoInfoG should be the filtered value
    collectMat(iProj,2) = collectMat(iProj,1)-FP.filoCount+FN.filoCount; % ground truth... 
    collectMat(iProj,3) = FP.filoCount; % false pos
    collectMat(iProj,4) = FN.filoCount;% false neg
    collectMat(iProj,5)=  collectMat(iProj,1) - collectMat(iProj,3); % True Positive: total number of detect - false positive 
    
    %NTot = sum(toKeep); 
    
%     nTrue(iProj,1) = NTot -FP.filoCount + FN.filoCount;
%     truePos = NTot - FP.filoCount; 
%     nPos(iProj,1) = FP.filoCount; 
%     nNeg(iProj,1) = FN.filoCount; 
%     percentPositive(iProj,1) = FP.filoCount/nTrue(iProj,1)*100; 
%     percentNegative(iProj,1) =FN.filoCount/nTrue(iProj,1)*100;

%     falseDis(iProj,1) = FP.filoCount/NTot; % false discovery rate : 
%     % ratio false positives to total detection 
%     detectSens(iProj,1) = (NTot-FP.filoCount)/nTrue(iProj,1); % proportion correctly IDed. 
%     prec(iProj,1) =   truePos/           % true positives/true positives + false positives
%     recall(iProj,1) = truePos/FN.filoCount+truePos; 
    



   if strcmpi(typeC,'branch')
       Mis = load([reconstructDir filesep toTest{iProj} filesep ...
           'Overlay' typeC filesep 'ValidationOverlays' filesep ...
           'Misconnection_' typeC '_Coords.mat']);
       collectMat(iProj,6) = Mis.filoCount;
       %         nmis(iProj,1) = Mis.filoCount;
       %         percentMis(iProj,1) = Mis.filoCount/nTrue(iProj,1).*100;
        
    end
    
    
end
% rowNames = toTest; 
%test = [falseDis detectSens nTrue nPos nNeg ]; 

% varNames{1,1} = ['False_Discovery_Rate_' typeC]; 
% varNames{1,2} = ['Detection Sensitivity_' typeC]; 
% varNames{1,3} = ['Calculated_Ground_Truth' typeC];
% 
% varNames{1,4} = ['Number_False_Positive_Filopodia' typeC]; 
% varNames{1,5} = ['Number_False_Negative_Filopodia' typeC]; 

% varNames{1,1} = ['Percent_False_Positive_' typeC]; 
% varNames{1,2} = ['Percent_False_Negative_' typeC]; 
% varNames{1,3} = ['Ground_Truth_' typeC 'Filopodia'];
% varNames{1,4} = ['Number_False_Positive_Filopodia' typeC]; 
% varNames{1,5} = ['Number_False_Negative_Filopodia' typeC]; 

% collect by group 
% test = cellfun(@(x) NTot(strcmpi( ))
% make a cell of values for each group 
 collectMatByGroup = cellfun(@(x)  collectMat(strcmpi(groupIDs,x),:),sampled,'uniformoutput',0); 
 
 rawValues{iType} = collectMatByGroup; 
 
 % precision : TP/TotalDetection (ie TP +FP) : Also known as positive
 % predictive power : want = to 1 
 %precisionByGroup = cellfun(@(x)  sum(x(:,5))./sum(x(:,1)),collectedMatByGroup,'uniformoutput',0);
 
 % false Discovery By Group : FP/NDetected (TP+FP) perfect = 0
 out{count} = cellfun(@(x) sum(x(:,3))./sum(x(:,1)),collectMatByGroup);
 calcName{count} = ['False Discovery Rate ' typeC];
 count = count+1;
 
 % false negative rate : FN/NTrue (TP+FN) perfect = 0
 out{count} = cellfun(@(x) sum(x(:,4))./sum(x(:,2)),collectMatByGroup);
 %
 calcName{count} = ['False Negative Rate ' typeC];
 count = count+1;
 
 
 if strcmpi(typeC,'branch')
     % branch misconnection rate :  Number Mis/Total TP - each row is
     % a group value, count is the calculation type 
     out{count} = cellfun(@(x) sum(x(:,6))./sum(x(:,5)),collectMatByGroup);
     calcName{count} = '% Misconnected Branches';
     count = count+1;
 end
 
 
 %recallByGroup = 

 % perform the calcs 
%  falseDis = cellfun(@(x) x(:,2)./x(:,1)    ;%  FP.filoCount/NTot
%  detectSens =  ; 
%  falseDisByGroup = cellfun(@(x) sum(x(:,3))./sum(x(:,1)),collectMatByGroup,'uniformoutput',0); 
%  detectSensByGroup = cellfun(@(x) sum(x(:,1)-sum(x(:,3))./sum(x(:,2))),collectMatByGroup,'uniformoutput',0);  % %(totalDetect-falsePositive)/true number

 % Recall : True Positive / (True Positive + FalseNegative) (ie ground
 % truth)
 %recallByGroup = cellfun(@(x)  sum(x(:,5))./(sum(:,5) + sum(:,4)),collectMatByGroup,'uniformoutput',0); % true positive/(true positive + false neg)-  
 
 % false negative rate : FN/(TP+FN) : perfect = 0 
 
 % false discovery rate : FP/(TP+FP) : FP/Total Detections  
 
% x = dataset({test,varNames{:}},'ObsNames',toTest); 
% export(x,'File',['PercentError' typeC '.csv']); 
clear collectMatByGroup
end % iType

% Make the HeatMap 

setAxis('on')
% compile the calculations (for a heat map) : rows are groups, columns calculation type
forHM = horzcat(out{:});

% calculate number of number of Cells, Frames, Filo  
% names = arrayfun(@(x)  ['Frames' num2str(size(collectMatByGroup{x},2)) 'Filo' sum(collectMatByGroup{x}(:,1))]



% Reformat so that control is last (on top of heat map)
idxControl = cellfun(@(x) ~isempty(regexp(x,'KDNo','ONCE')),sampled); 

frameNum = cellfun(@(x) size(x,1) , rawValues{1}); 
sampled = arrayfun(@(x) [sampled{x} ' Frames ' num2str(frameNum(x))],1:length(frameNum),'uniformoutput',0); 

sampled = [sampled(~idxControl) sampled(idxControl)]; 

sampled = cellfun(@(x) strrep(x,'KD',''),sampled,'uniformoutput',0); 

sampled = cellfun(@(x) strrep(x,'No','Control'),sampled,'uniformoutput',0); 
sampled = cellfun(@(x) strrep(x,'Dock','Dock7'),sampled,'uniformoutput',0); 
sampled = cellfun(@(x) strrep(x,'CDC','CDC42'),sampled,'uniformoutput',0); 

% add frame number 
%frameNum = cellfun(@(x) size(

forHM = [ forHM(~idxControl,:);forHM(idxControl,:)];

%cmap = brewermap(128,'OrRd');
imagesc(forHM);
hold on 
%colormap(cmap);
colormap(parula)
%  HeatMap(forHM,'RowLabels',sampled,'ColumnLabels',calcName, 'colormap','redbluecmap','ColumnLabelsRotate', ...
%     45);
set(gca,'YTick',1:numel(sampled));

set(gca,'YTickLabel',sampled);

set(gca,'XTick',1:numel(calcName));
set(gca,'XTickLabel',calcName);
set(gca,'XTickLabelRotation',45); 
colorbar

% format just that Control is on top 

% arrayfun(@(x,y) text(x,y, ['Frames: ' num2str(rawValues{y}{x},2) ' Filo: ' ...
%     numstr(rawValues{y}{x}(:,1))],'Color','w'),1:numel(calcName),1:numel(sampled));   

% for now just do a quick fix for N 
NFilo = arrayfun(@(y) cellfun(@(x) sum(x(:,1)),rawValues{y}),1:3,'uniformoutput',0);  
NFilo = cellfun(@(x) [x(~idxControl); x(idxControl)],NFilo,'uniformoutput',0); 
pos = [1,3,6]; 

%arrayfun(@(x,y) text(pos(x),y,['N Filo:' num2str(NFilo{x}(y))],'color','w'),1:3,1:numel(sampled)); 

for iType = 1:3 
    for iGroup= 1:numel(sampled)
    text(pos(iType)-0.3,iGroup,['N Filo ' num2str(NFilo{iType}(iGroup))],'color','r','FontSize',10); 
    end     
end 

saveas(gcf,[outPutDirectory filesep 'ErrorHeatMap.fig']);
saveas(gcf,[outPutDirectory filesep 'ErrorHeatMap.eps'],'psc2'); 
saveas(gcf,[outPutDirectory filesep 'ErrorHeatMap.png']); 


calcNameForTable = cellfun(@(x) strrep(x,' ','_'),calcName,'uniformoutput',0); 
% t = array2table(forHM ,'VariableNames',calcNameForTable,'RowNames',sampled');
% writetable(t,[outPutDirectory filesep 'ErrorTable.csv']); 
d = mat2dataset(forHM,'VarNames',calcNameForTable,'ObsNames',sampled); 
export(d,'file',[outPutDirectory filesep 'ErrorTable.csv']); 

save([outPutDirectory filesep 'RawValues'],'rawValues'); 

% False_PositiveEmbedded_Percent = percentPositive.*100; 
% False_NegativeEmbedded_Percent = percentNegative.*100;  

% t  = table([False_PositiveEmbedded_Percent,False_NegativeEmbedded_Percent],'RowNames',rowNames); 
% writetable(t,'PercentErrorEmbedded.csv'); 
end 
function [groupID] = getGroupIDs(toTest)
%name = upDirectory(toTest,1); 

% get ID from name (Note should make more generic) 
name = strrep(toTest,'_','/');
[~,groupID] = upDirectory(name,4); 
end 
