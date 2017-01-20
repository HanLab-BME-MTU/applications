function [matPosHist,matNegHist] = makeDescriptorMIP(matPosMIP,matNegMIP)
%Makes Histograms given positive and negative MIP stacks

%Initialize the histograms (dim: 256 bins x 8 channels)
matPosHist = zeros(256,8);
matNegHist = zeros(256,8);

%Loop over the 8 alpha channels
for alpha = 1:8
    %Linearize the matrix
    curvecpos = reshape(matPosMIP(:,:,alpha),1,[]);
    curvecneg = reshape(matNegMIP(:,:,alpha),1,[]);

    %Calculate the histograms
    temp_histNeg = histc(curvecneg,0:255);
    temp_histPos = histc(curvecpos,0:255);
    
    %Force the counts of 0 to zero
    temp_histNeg(1) = 0;
    temp_histPos(1) = 0;
   
    %Store values
    matPosHist(:,alpha) = temp_histPos;
    matNegHist(:,alpha) = temp_histNeg;
end
end