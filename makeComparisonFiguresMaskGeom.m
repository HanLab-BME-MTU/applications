% function makeComparisonFiguresMaskGeom(conStructs,conNames,isGood)
% 
% nRes = numel(conStructs);
% 
% sToHist = { 'mcMeanPerMov'
%             'mcMedPerMov'
%             'mcVarPerMov'
%             'mapcMeanPerMov'
%             'mapcMedPerMov'
%             'mapcVarPerMov'};
% binsHist = {'histBinsPerMov'
%             'histBinsPerMov'
%             'histBinsPerMov'
%             'mapcHistBinsPerMov'        
%             'mapcHistBinsPerMov'        
%             'mapcHistBinsPerMov'};
%     
% conCols = jet(nRes);        
%         
% for j = 1:numel(sToHist)
%     
%     figure
%     hold on
%     for k = 1:nRes
%         n(:,k) = histc(horzcat(conStructs(k).(sToHist{j}){isGood{k}}),conStructs(k).(binsHist{j}){1});
%         n(:,k) = n(:,k) ./ sum(n(:,k));
%     end
%     bar(n,repmat(conStructs(1).(binsHist{j}){1})
%     legend(conNames)
%     xlabel(sToHist{j})
%     ylabel('Normalized Histogram')
% end
%     
%     
    
%NOte - two outliers in WT with incorrect pixel size (#s 1 and 2) are
%excluded from analyiss

outDir = 'L:\nih\Post Processing Myo Inhib and WT\Mask Geometry\Comparisons';

nBins = 10;

figure
bins = linspace(6e-4,9e-4,nBins);
nWT = histc(horzcat(wt.mapcMedPerMov{3:end}),bins);
nWT = nWT ./ sum(nWT);
nBleb = histc(horzcat(bleb.mapcMedPerMov{:}),bins);
nBleb = nBleb ./ sum(nBleb);

bar([bins' bins'],[nWT' nBleb'],1)
    
xlabel('Median of Maximum Absolute Curvature Component, 1/nm')
ylabel('Normalized Count')
legend('Untreated','Blebbistatin')

saveThatShit('Median Max Abs Curv comparison wt and bleb bar graph',outDir);

plot(bins',nWT', bins', nBleb','LineWidth',3)
    
xlabel('Median of Maximum Absolute Curvature Component, 1/nm')
ylabel('Normalized Count')
legend('Untreated','Blebbistatin')

saveThatShit('Median Max Abs Curv comparison wt and bleb line plot',outDir);


figure
bins = linspace(-3e-4,-1e-4,nBins);
nWT = histc(horzcat(wt.mcMedPerMov{3:end}),bins);
nWT = nWT ./ sum(nWT);
nBleb = histc(horzcat(bleb.mcMedPerMov{:}),bins);
nBleb = nBleb ./ sum(nBleb);

bar([-bins' -bins'],[nWT' nBleb'],1)
    
xlabel('Median of Mean Curvature, 1/nm')
ylabel('Normalized Count')
legend('Untreated','Blebbistatin')

saveThatShit('Median Mean Curv comparison wt and bleb bar plot',outDir);

plot(-bins',nWT', -bins', nBleb','LineWidth',3)

xlabel('Median of Mean Curvature, 1/nm')
ylabel('Normalized Count')
legend('Untreated','Blebbistatin')

saveThatShit('Median Mean Curv comparison wt and bleb line plot',outDir);




