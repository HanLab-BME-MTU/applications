function [res] = determineDetectedObjectSizes(exp)
% determine size of clathrin pits in all images for a given condition
%
% SYNOPSIS:
% [pitResults]=determineAllPitsizes(exp);
%
% INPUT     :   exp   = field containign all exp for the condition
%
% OUTPUT    :   res = structure containing for every movie
%                   .sigma : sigma of Gaussian fit
%                   .A     : amplitude
%                   .c     : background
%
%
% REMARKS   :
%
% created with MATLAB ver.: 7.9.0.529 (R2009b) on MacOS X 10.6
%
% Francois Aguet, 01/24/2010

nExp = length(exp);
data(1:nExp) = struct('image', [], 'detection', []);

for i = 1:nExp
    
    % load list of frames
    frameList = dir([exp(i).source '*.tif']);
    frameList = {frameList.name};
    
    if isempty(frameList)
        framePath = [uigetdir(exp(i).source, 'Select directory containing the movie frames:') filesep];
        frameList = dir([framePath '*.tif']);
        frameList = {frameList.name};
        if isempty(frameList)
            error(['No frames found in ' exp(i).source]);
        end
    end
    data(i).image = double(imread([exp(i).source frameList{1}]));
    
    % load 'detection.mat'
    detectionFile = load([exp(i).source 'DetectionStructures' filesep 'detection.mat']);
    if isfield(detectionFile,'detection')
        data(i).y = detectionFile.detection(1).xCoord(:,1);
        data(i).x = detectionFile.detection(1).yCoord(:,1);
    elseif isfield(detectionFile,'cdet')
        data(i).y = detectionFile.cdet(1).xCoord(:,1);
        data(i).x = detectionFile.cdet(1).yCoord(:,1);
    else
        error('no detection exp file of specified format found');
    end
end

res(1:nExp) = struct('sigma', [], 'A', [], 'c', []);
for i=1:nExp
    fprintf('Processing movie no. %d\n',i);
    res(i) = determinePitsizeFromImage_2d(data(i).image, data(i).x, data(i).y);
end

mu = mean(res.sigma);
d = 2*mu/150;
xa = 0:d:2*mu;

[nv, xv] = hist(res.sigma, xa);
figure; bar(xv, nv);
axis([0 xa(end) 0 max(nv)]);
set(gca, 'FontName', 'Helvetica', 'FontSize', 14, 'LineWidth', 1.5);
title(['Sigma: mean = ' num2str(mu, '%3.3f') ', std = ' num2str(std(res.sigma), '%3.3f')]);

end % of function


function [results] = determinePitsizeFromImage_2d(image, x, y)
% determine size of clathrin pits in image (in pixels)
%
% SYNOPSIS:
% [results]=determinePitsizeFromImage(image, xypos);
%
% INPUT     :   image   = tiff image containing pits
%               xypos   = nx2 vector containing x,y positions of pits,
%                         which can be extracted from detection results
%
% OUTPUT    :   results = vector containing the results of the Gaussian fit
%                         to the pits, with
%                         [Background Amplitude Gausswidth];

ru = 5;
[nx, ny] = size(image);
N = length(x);

sigma = zeros(1,N);
A = zeros(1,N);
c = zeros(1,N);
valid = zeros(1,N);

for k = 1:N
    
    x0 = round(x(k));
    y0 = round(y(k));
    
    % ignore objects at image boundaries
    if (x0-ru>=1 && y0-ru>=1 && x0+ru<=nx && y0+ru<=ny)
        window = image(x0-ru:x0+ru, y0-ru:y0+ru);

        imax = double(max(window(:)));
        imean = double(mean(window(:)));
        
        % Gaussian parameter vector: [x y A sigma background]
        
        % sigma fixed, other parameters free
        [estimates1] = fitGaussian2D(window, [ru+1 ru+1 imax-imean 1.5 imean], 'xyAc');
        
        % background and positions fixed, estimate A and sigma
        [estimates2] = fitGaussian2D(window, estimates1, 'As');
        
        sigma(k) = estimates2(4);
        A(k) = estimates2(3);
        c(k) = estimates2(5);
        valid(k) = 1;
        if (A(k) < 0 || sigma(k) < 0)
            valid(k) = 0;
        end;
    end
end
sigma(~valid) = [];
A(~valid) = [];
c(~valid) = [];
results.sigma = sigma;
results.A = A;
results.c = c;
end