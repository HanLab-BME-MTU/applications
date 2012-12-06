function plotDataSetsResult(dataSet,label)

if ~isfield(dataSet,{'meanValue','CI'})
    error('Non-recongnizable Structure')
end

nSet      = numel(dataSet);
condition = fieldnames(dataSet(1).CI);
nCond     = numel(condition);
variable  = fieldnames( dataSet(1).CI.(condition{1}) );
nVar      = numel(variable);

for iVar = 1:nVar
    
    for iSet = 1:nSet
        
        for iCond = 1:nCond
            
            meanDisp(iSet,iCond,iVar) = dataSet(iSet).meanValue.(condition{iCond}).(variable{iVar});
            UpCI(iSet,iCond,iVar)     = dataSet(iSet).CI.(condition{iCond}).(variable{iVar})(2);
            LwCI(iSet,iCond,iVar)     = dataSet(iSet).CI.(condition{iCond}).(variable{iVar})(1);
            UpCI(iSet,iCond,iVar)     = UpCI(iSet,iCond,iVar) - meanDisp(iSet,iCond,iVar);
            LwCI(iSet,iCond,iVar)     = meanDisp(iSet,iCond,iVar) - LwCI(iSet,iCond,iVar);
            
        end
        
    end
    
    figure
    barplot2(meanDisp(:,:,iVar), UpCI(:,:,iVar), LwCI(:,:,iVar),'ErrorBarPosition', 'both',...
      'BarWidth', 0.8, 'GroupDistance', 2,'Xlabels',label.dataSet,'YLabel', label.variable{iVar});
    h1 = legend(label.condition);
    set(h1,'EdgeColor','w')
   
end