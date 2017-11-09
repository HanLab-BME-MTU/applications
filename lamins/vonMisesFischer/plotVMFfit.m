function [ out ] = plotVMFfit( x )
%plotVMFfit plot vonMisesFischer distribution using fitting parameters

N = 36;
theta = -pi/2:pi/N:pi/2-pi/N;

vmf = vonMisesFischer2d(theta);
hold on;
h1 = plot(theta,x(2:2:end-2)*vmf(x(1:2:end-2)',x(end-1)) + x(end));
h2 = plot(theta,(repmat(x(2:2:end-2)',1,N).*vmf(x(1:2:end-2)',x(end-1)) + x(end))');

out = [h1; h2];


end

