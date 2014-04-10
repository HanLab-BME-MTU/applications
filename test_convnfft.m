
clc
clear

fprintf( '\n\nComparing the performance of convolution in spatial vs frequence domain ... \n\n' );
A = zeros(512,512,40);

sigma = 6;
w = 5*sigma;
B = ones(w,w,w);
B = B / numel(B);

tic
C = convnfft(A,B,'same');
convnfft_Time = toc

tic
C = convn(A,B,'same');
convn_Time = toc

freqSpeedup = convn_Time / convnfft_Time