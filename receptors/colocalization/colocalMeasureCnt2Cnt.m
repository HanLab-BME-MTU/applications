function [r, randR] = colocalMeasureCnt2Cnt(imageCntA,imageCntB,maskingFile)
% COLOCALMEASURECNT2CNT measures colocalization for two channels where both are continuous
%

    % Read in continuum image
    ICntA = double(imageCntA); 
    ICntB = double(imageCntB); 
    
    %Correct images for non-uniform background
    compValue = mean(ICntA(maskingFile~=0));
    compMask = compValue*ones(256,256);
    nImage = filterGauss2D(ICntA,10);
    ICntA = ICntA-nImage;
    ICntA = ICntA+ compMask;
    
    compValue = mean(ICntB(maskingFile~=0));
    compMask = compValue*ones(256,256);
    nImage = filterGauss2D(ICntB,10);
    ICntB = ICntB-nImage;
    ICntB = ICntB+ compMask;
    
    
    r = corr2(ICntA,ICntB);
    [~,ImageRand] = RandomizeBlocks(ICntA,maskingFile,[16,16]);
    randR = corr2(ImageRand,ICntB);
end