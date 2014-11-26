function [matDist] = analyzePlasticCell(matPosHist,matNegHist,all_size,vecFragmentVideoNumLabel,SELECTION_NUM)
%normalize features respect to image size
for jj = 1:length(all_size(:,1))
    matPosHist(jj,:,:) = matPosHist(jj,:,:)./all_size(jj,1);
    matNegHist(jj,:,:) = matNegHist(jj,:,:)./all_size(jj,1);
end


matCellPosHist = matPosHist(vecFragmentVideoNumLabel == SELECTION_NUM,:,:);
matCellNegHist = matNegHist(vecFragmentVideoNumLabel == SELECTION_NUM,:,:);

matPosDist = zeros(size(matCellPosHist,1),size(matCellPosHist,1));
matNegDist = zeros(size(matPosDist));
for i = 1:size(matCellPosHist,1)
    for j = 1:size(matCellPosHist,1)
        matPosDist(i,j) = sum(abs(matCellPosHist(i,:,:) - matCellPosHist(j,:,:)),3);
        matNegDist(i,j) = sum(abs(matCellNegHist(i,:,:) - matCellNegHist(j,:,:)),3);
    end
end
matDist = matPosDist + matNegDist;
end