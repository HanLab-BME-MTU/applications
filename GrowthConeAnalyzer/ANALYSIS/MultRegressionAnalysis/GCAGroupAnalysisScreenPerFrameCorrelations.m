function [ toPlot ] = GCAGroupAnalysisScreenPerFrameCorrelations(toPlot)
% GCAGroupAnalysisScreenPerFrameCorrelations


% collect the fields to collect 
outDir = pwd; 
measNames = fieldnames(toPlot);

measNames = measNames(~(strcmpi(measNames,'functional') | strcmpi(measNames,'info') | strcmpi(measNames,'protrusionAnalysis_persTime_greaterThan75th_Perc') ... 
    | strcmpi(measNames,'retractionAnalysis_persTime_greaterThan75th_Perc'))); 



elongVelCell = toPlot.functional.neuriteElongationVel;
% first just look at pooled metrics

nGroups = 1;

for iGroup = 1:nGroups
    
    elongVelC = elongVelCell{iGroup} ;
    % crop each by 119
    elongVelCCrop = cellfun(@(x) x(1:119)',elongVelC,'uniformoutput',0);
    
    outDirG = [outDir filesep toPlot.info.names{iGroup} 'PerTrajectoryCorrMatrix']; 
    if ~isdir(outDirG)
        mkdir(outDirG); 
    end 
    
    nProjs = size(toPlot.info.projList{iGroup},1); 
    %     elongVelAll = vertcat(elongVelCCrop{:});
    for iProj = 1:nProjs
        
        projID = gcaGetNeuriteID([toPlot.info.projList{iGroup}{iProj,1} filesep 'GrowthConeAnalyzer' ]);
        
        projID = strrep(projID,' ','_');
        % reformat the dataMat per cell such that each meas a column and each row a frame.
        % as currently the data is stored such that each meas is a
        % structure field, each row is a frame, and each column is a
        % project...
        dataPerMeas   = cellfun(@(x) [ toPlot.(x).dataMatPerFrame{iGroup}(:,iProj)],measNames,'uniformoutput',0);
        dataPerMeas = dataPerMeas'; 
        dataPerMeasFinal = [ horzcat(dataPerMeas{:}) elongVelCCrop{iProj}]; 
        
        % unfortunately need to remove NaNs. 
        
        idxRemove = sum(isnan(dataPerMeasFinal),2); 
        dataPerMeasFinal = dataPerMeasFinal(~idxRemove,:); 
        
        toPlot.corr.dataPerMeas{iGroup}{iProj} = dataPerMeasFinal;
        [rAll,pAll] = corr(dataPerMeasFinal,'type','Spearman');
        toPlot.corr.rAll{iGroup}{iProj} = rAll;
        toPlot.corr.pAll{iGroup}{iProj} = pAll;
       
        
        
        
        
        % screen pAll 
        pNEvsS = pAll(end,1:end-1); 
        measSig = measNames(pNEvsS<0.05); 
        toPlot.corr.measSig{iGroup}{iProj} = measSig; 

        idxToPlot = find(pNEvsS<0.05); 
        
        if ~isempty(measSig)
            
          for i = 1:length(idxToPlot)
              cDir = [outDirG filesep measSig{i}];
              if ~isdir(cDir)
                  mkdir(cDir); 
              end 
              
              setAxis('off')
              scatter(dataPerMeasFinal(:,idxToPlot(i)),dataPerMeasFinal(:,end) ,50,'k'); 
              title({['Spearman Correlation Coef = ' num2str(rAll(end,idxToPlot(i)),3)] ; ... 
                  ['PValue = ' num2str(pAll(end,idxToPlot(i)),3)]}); 
              xlabel(measSig{i}); 
              ylabel({'Instantaneous Neurite Elongation Velocity'; 'um/min'}); 
              saveas(gcf,[cDir filesep projID '.fig']); 
              saveas(gcf,[cDir filesep projID '.eps'],'psc2');
              saveas(gcf,[cDir filesep projID '.png'])
              close gcf
          end 
          
        end 
        
        
    end
    
    
    % collect all info per group 
    dataPerMeasAll = vertcat(toPlot.corr.dataPerMeas{iGroup}{:}); 
    [rAll,pAll] = corr(dataPerMeasAll,'type','Spearman'); 
    
   
        
        % screen pAll 
        pNEvsS = pAll(end,1:end-1); 
        measSig = measNames(pNEvsS<0.05); 
        toPlot.corr.measSig{iGroup}{iProj} = measSig; 
        idxToPlot = find(pNEvsS<0.05); 
         save([outDirG filesep 'corrPerFrameSum.mat'],'dataPerMeasAll','rAll','pAll','measSig'); 
         save([outDirG filesep 'toPlotWithCorr.mat'],'toPlot'); 
        if ~isempty(measSig)
            
          for i = 1:length(idxToPlot)
              cDir = [outDirG filesep measSig{i}];
              if ~isdir(cDir)
                  mkdir(cDir); 
              end 
              
              setAxis('off')
              scatter(dataPerMeasAll(:,idxToPlot(i)),dataPerMeasAll(:,end) ,50,'k'); 
              title({['Spearman Correlation Coef = ' num2str(rAll(end,idxToPlot(i)),3)] ; ... 
                  ['PValue = ' num2str(pAll(end,idxToPlot(i)),3)]}); 
              xlabel(measSig{i}); 
              ylabel({'Instantaneous Neurite Elongation Velocity'; 'um/min'}); 
              saveas(gcf,[cDir filesep  'All.fig']); 
              saveas(gcf,[cDir filesep  'All.eps'],'psc2');
              saveas(gcf,[cDir filesep  'All.png'])
              close gcf
          end 
          
        end 
        
end


end

%     for iMeas = 1:numel(measNames)
%         
%         % Test Each Plot 
%         
%         dataMatC = toPlot.(measNames{iMeas}).dataMatPerFrame{iGroup};
%         
%         [r,p ] = arrayfun(@(x) corr(elongVelCCrop{x},dataMatC(:,x),'type','Spearman'),1:20); 
%          
%         toPlot.(measNames{iMeas}).rp
%         
%         % make a folder for the 
%         
%        %[r,p] = cellfun(@(x,y) corr([x ,y],'name','Spearman'),elongVelCCrop,dataMatC); 
%         
%        h = arrayfun(@(i) scatter(elongVelCCrop{i},dataMatC{i},50,'k','filled'),1:20,'uniformoutput',0); 
%         
%     end
    



