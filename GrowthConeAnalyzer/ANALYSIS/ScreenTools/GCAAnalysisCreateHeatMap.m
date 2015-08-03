function [ dataMat ] = GCAAnalysisCreateHeatMap( toPlot ,displayRange)
%
if nargin<2 
    displayRange = []; 
end 
% collect mean values for the per neurite data per group 

fields = fieldnames(toPlot); 
% take out info 
fields = fields(cellfun(@(x) ~strcmpi(x,'info'),fields)); 

for iParam = 1:numel(fields)
    % for now do this per cell : each column of the toPlot.param{iGroup} 
    % contains a row of observations for that neurite
    % get the mean of the parameter for each group 
    
      meansAllGroups = cellfun(@(x)  nanmean(nanmedian(x,1)),toPlot.(fields{iParam})); 
      
      % for now just get the percent change 
      perChange{iParam,1} = (meansAllGroups-meansAllGroups(1))./meansAllGroups(1)*100; 
     
    % means
    

end

% add two extra fields for the median of the 75th percentile of protrusion
% persistence and retraction persistence. 
 

fieldsFinal = cell(numel(fields)+3,1); 
fieldsFinal(1:numel(fields))= fields; 
fieldsFinal{end-1} = ''; 
fieldsFinal{end-2} = '';
fieldsFinal{end} = 'Dummy'; 

dataMat = vertcat(perChange{:}); % make the  data Mat such that each row is a 
 
 
 
 yLabels = fieldsFinal; 
 xLabels = toPlot.info.names; 
%  displayRange  = repmat(displayRange,1,(numel(xLabels))); 
%  buffer = zeros(2,numel(xLabels)); 
%  dataMat = [dataMat;buffer;displayRange]; 
 %dataMat = dataMat.*100;
 colorMap = HeatMap(dataMat(:,2:end),'RowLabels',fieldsFinal,'ColumnLabels',toPlot.info.names(2:end),... 
    'colormap','redbluecmap','ColumnLabelsRotate',45); 
save('colorMap','colorMap'); 
 %colorbar
%gca_heatmap(dataMat',yLabels,xLabels,[]); 
