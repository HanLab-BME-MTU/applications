function  [] = displayColocalizationResults(conditions,resultType)
%General Molecule Enrichment
%Load colocalization Data
% for k = 1:length(conditions)
%     load(strcat(conditions{k},'/colorStackMovieData.mat'))
%     p = MD.getProcess(3).getParameters();
%     p.ChannelObs = 1;
%     MD.getProcess(3).setParameters(p);
%     MD.getProcess(3).run;
% end

if resultType == 1
for m = 1:(length(conditions))
load(strcat(conditions{m},'/colocalInfoCT.mat'),'randRatioAve')
% load(strcat(conditions{m},'/colocalInfo.mat'),'randRatioAve')
molecules(:,m) = randRatioAve(1:30,1);
% molecules(:,m+1) = randRatioAve(:,1);

% molecules2(:,m) = ratioAve(:,2);
% molecules2(:,m+1) = randRatioAve(:,2);

end
figure; boxplot(molecules,'notch','on','labels',conditions);
% figure; boxplot(molecules2,'notch','on');
elseif resultType == 2
    for m = 1:length(conditions)
        load(strcat(conditions{m},'/colocalInfo.mat'),'cellIntensity')
        molecules(:,m) = cellIntensity(:,1);


    end 
    figure; boxplot(molecules,'notch','on','labels',conditions);
    title('Untreated')
end
%Specific Enrichment in Actin Regimes

% molecules(:,1) = cellCtHighU1(:,1);
% molecules(:,2) = cellCtLowU1(:,1);
% molecules(:,3) = cellCtHighT1(:,1);
% molecules(:,4) = cellCtLowT1(:,1);
% molecules(:,5) = cellCtHighG1(:,1);
% molecules(:,6) = cellCtLowG1(:,1);
% molecules(:,7) = cellCtHighM1(:,1);
% molecules(:,8) = cellCtLowM1(:,1);
% figure; boxplot(molecules,'notch','on');
% 
% ax = gca;
% set(ax,'XTickLabel',{'NT High','NT Low','+TSP High','+TSP Low','IgG High','IgG Low','IgM High','IgM Low'})
% 
% molecules(:,1) = cellCtHighU(:,2);
% molecules(:,2) = cellCtLowU(:,2);
% molecules(:,3) = cellCtHighT(:,2);
% molecules(:,4) = cellCtLowT(:,2);
% molecules(:,5) = cellCtHighG(:,2);
% molecules(:,6) = cellCtLowG(:,2);
% molecules(:,7) = cellCtHighM(:,2);
% molecules(:,8) = cellCtLowM(:,2);
% figure; boxplot(molecules,'notch','on');
% 
% ax = gca;
% set(ax,'XTickLabel',{'NT High','NT Low','+TSP High','+TSP Low','IgG High','IgG Low','IgM High','IgM Low'})
% %-----------------------------------------------------------------------------------------
% molecules(:,1) = ratioAveU(:,1);
% molecules(:,2) = randRatioAveU(:,1);
% molecules(:,3) = ratioAveT(:,1);
% molecules(:,4) = randRatioAveT(:,1);
% molecules(:,5) = ratioAveG(:,1);
% molecules(:,6) = randRatioAveG(:,1);
% molecules(:,7) = ratioAveM(:,1);
% molecules(:,8) = randRatioAveM(:,1);
% figure; boxplot(molecules,'notch','on');
% ax = gca;
% set(ax,'XTickLabel',{'NT','NT Rand','+TSP','+TSP Rand','IgG','IgG Rand','IgM','IgM Rand'})
% 
% molecules(:,1) = ratioAveU(:,2);
% molecules(:,2) = randRatioAveU(:,2);
% molecules(:,3) = ratioAveT(:,2);
% molecules(:,4) = randRatioAveT(:,2);
% molecules(:,5) = ratioAveG(:,2);
% molecules(:,6) = randRatioAveG(:,2);
% molecules(:,7) = ratioAveM(:,2);
% molecules(:,8) = randRatioAveM(:,2);
% figure; boxplot(molecules,'notch','on');
% ax = gca;
% set(ax,'XTickLabel',{'NT','NT Rand','+TSP','+TSP Rand','IgG','IgG Rand','IgM','IgM Rand'})
% 
% molecules(:,1) = ratioAveU(:,1);
% molecules(:,2) = ratioAveU1(:,1);
% molecules(:,3) = randRatioAveU(:,1);
% molecules(:,4) = randRatioAveU1(:,1);
