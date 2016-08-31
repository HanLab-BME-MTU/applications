function samplingSimulation(interval,noise_ratio)

close all;
x = 0:.1:30;
y = normpdf(x,15,1);
figure;
plot(x,y);
axis([0 30 0 0.5]);
title('Single Point Spread Function');

% Interval value should be integral multiple of 0.01 and no larger than 30
index = 1 + floor(interval*rand/.1);
step = interval/0.1;
sampleSize = floor((length(x)-index)/step);
xSampling = zeros(1,sampleSize);
ySampling = zeros(1,sampleSize);
yNoise = zeros(1,sampleSize);
for n = 1:sampleSize
    xSampling(n) = x(index);
    ySampling(n) = y(index);
    
    % Add noise, suggestted 10-30% of Gaussian distribution amplitude
    yNoise(n) = y(index) + noise_ratio*max(y)*rand(1);
    index = index + step;
end
figure;
plot(xSampling,ySampling);
axis([0 30 0 0.5]);
title('PSF after sampling');

% Plot distribution with noise
figure;
plot(xSampling,yNoise);
axis([0 30 0 0.5]);
title('Sampled PSF with noise');

% Plot distribution after gaussian filting
filter = fspecial('gaussian',6,.5);
yfilt = imfilter(yNoise,filter);
figure;
plot(xSampling,yfilt);
axis([0 30 0 0.5]);
title('Sampled PSF after Gaussian Filtering');