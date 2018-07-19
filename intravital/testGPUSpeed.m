
clc
clear all
close all

A = ones(512,512,60);
B = ones(512,512,60);

fprintf(1, '\nRunning ion CPU...\n');
tic
C = A + A .* B - 2 * B;
toc

fprintf(1, '\nRunning ion GPU...\n');
tic
AG = gpuArray(A);
BG = gpuArray(B);
CG = AG + AG .* BG - 2 * BG;
C = gather(CG);
toc