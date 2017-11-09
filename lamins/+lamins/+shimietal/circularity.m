D = dir('/project/biophysics/jaqaman_lab/lamins/2015/20150602/MEF*');
circularity = arrayfun(@(s) nanmean(vertcat(s.nucleusCircularity.maskCircularity)),stats)';
circularityTable = cell2table([{D.name}' num2cell(circularity)],'VariableNames',{'Set','Circularity'});
circularityTable.NinetyNine_FaceArea = arrayfun(@(s) prctile(s.area(1).all.data,99),stats)';