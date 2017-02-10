function [ out ] = analyzeExamples( I )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

if(nargin < 1)
    D = dir('*.png');
    I = {D.name};
end

if(iscellstr(I))
    parcellfun_progress(@intersections.analyzeExamples,I,'Unif',false);
    return;
end

if(ischar(I))
    filename = strrep(I,'.png','');
    I = imread(I);
else
    filename = [];
end

if(~isempty(filename))
    mkdir('nlms_vs_nms');
    mkdir('nlms_vs_I');
    mkdir('nlms3_m8_merge_vs_I');
    mkdir('nlms_gated_vs_I');
end

F = OrientationSpaceFilter.constructByRadialOrder(1/2/pi./2,1,8);
R = F*I;
R3 = R.getResponseAtOrderFT(3);
maxima = R.getRidgeOrientationLocalMaxima;
nlms3 = R3.nonLocalMaximaSuppression;
nlms3_m8 = R3.nonLocalMaximaSuppressionPrecise(maxima);
nms = nlms3(:,:,1);
nlms_mip = nanmax(real(nlms3_m8),[],3);
nlms3_m8_mip = nanmax(real(cat(3,nlms3,nlms3_m8)),[],3);

if(~isempty(filename))
    imwrite(imfuse(nms,nlms_mip),['nlms_vs_nms/' filename '_nlms_vs_nms.png']);
    imwrite(imfuse(nlms3_m8_mip,I),['nlms3_m8_merge_vs_I/' filename '_nlms_vs_I.png']);
    imwrite(imfuse(nlms_mip,I),['nlms_vs_I/' filename '_nlms_vs_I.png']);
    imwrite(imfuse(nlms_mip.*double(I),I),['nlms_gated_vs_I/' filename 'nlms_gated_vs_I.png']);
end


end


