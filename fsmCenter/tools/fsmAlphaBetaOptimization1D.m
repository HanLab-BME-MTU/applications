function [noiseParameter, actualP, Lmax_ratio] = fsmAlphaBetaOptimization1D(imageNameList, noiseParameter, confidenceP, bitDepth)
% Function alphaBetaOptimization1D optimizes noise parameters alpha & beta based on their given initial values
% such that the speckle detection result best matches user-specified confidence probability. The optimization
% search is performed in a 1D fashion, namely first with respect to beta, then with respect to alpha.
%
% SYNOPSIS   [noiseParameter] = alphaBetaOptimization1D(optimImageStack, optimImageStackLen, noiseParameter)
%
% INPUT        imageNameList      :     self-explanatory.
%             noiseParameter      :     original noise parameters in form [k sigmaDN PoissonNoise I0] 
%                                       where k = quantile / GaussianRatio;
%                confidenceP      :     user specified confidence probability. The objective of the optimization process
%                                       is to make the Percentage of ROI Speckles as close to this value as possible. 
% 
% ANNOTATION   Percentage of ROI Speckles is the percentage of accepted speckles within ROI (region of interest)
%              with respect to total accepted speckles.
%                              
% OUTPUT      noiseParameter      :     optimized noise parameters among which
%                                       sigmaD (alpha) and PossionNoise (beta) have been optimized. 
%
% DEPENDENCES     fsmPrepBkgEstimDelauNoEnh
                 

% Ge Yang, Aprile 8, 2004

% Step 1: --------------------------------------------------------
% Setting ROI through user input. Compute Gauss filtered image
% stack, local_max stack, local_min stack and cands stack

if (isempty(imageNameList)) % in case the name list is not provided
    [firstImageName, pathname] = uigetfile('.tif', 'Please choose the first image from the sequence for optimization.');
    firstImageName = strcat(pathname, firstImageName);
    imageNameList = getFileStackNames(firstImageName);
    optimImageStackLen = 3;  % use 3 frames by default
else
    firstImageName = imageNameList(1);
    optimImageStackLen = length(imageNameList);
end

info = imfinfo(char(firstImageName), 'tiff');
imHeight = info.Height;
imWidth = info.Width;

[img_temp, map_temp] = imread(char(firstImageName), 'tiff');
figure;
imshow(img_temp, map_temp);
bw = roipoly;
% imshow(bw);
% pause;
close;

% Gauss filtering, find localmax and localmin
imGStack = zeros(imHeight, imWidth, optimImageStackLen);   % Gauss filtered image stack
imMinStack = zeros(imHeight, imWidth, optimImageStackLen); % local minimum stack
imMaxStack = zeros(imHeight, imWidth, optimImageStackLen); % local maximum stack

h = waitbar(0, 'Computing local maxima & minima over the image stack, Please wait...');
cStack = struct('cands',[]); % Stack of cands
sigmaG = 1.0;

for i = 1 : optimImageStackLen
    img = double(imread(char(imageNameList(i)), 'tiff')) / (2 ^ bitDepth-1);  % Normalization
    imGStack(:, :, i) = gauss2d(img, sigmaG);
    imMinStack(:, :, i) = locmin2d(imGStack(:, :, i), [3,3]);
    imMaxStack(:, :, i) = locmax2d(imGStack(:, :, i), [5,5]);
    cStack(i).cands = fsmPrepBkgEstimDelauNoEnh(size(imGStack(:, :, i)), imMaxStack(:, :, i), imMinStack(:, :, i)); % Finds 3 loc min around each loc max
    waitbar(i/optimImageStackLen, h)
end
close(h);

h = msgbox('Now starting optimization search');
pause(2);
close(h);

total_area = imHeight * imWidth;
roi_area = sum(sum(bw))
non_roi_area = total_area - roi_area;

% fprintf('Percentage of roi area is %f\n', roi_area/total_area *100);
% pause;
% pause;

% Step 2: --------------------------------------------------------
% Setting ROI through user input. Compute Gauss filtered image
% stack, local_max stack, local_min stack and cands stack

% Step 2.1: Set grid search parameters.
beta_num = 10;            % number of beta search steps
beta_num_half = 0.5 * beta_num;
beta_search_step = 0.2;   % increase/decrease using stepsize of 20 percent

alpha_num = 4;            % number of alpha search steps
alpha_num_half = 0.5 * alpha_num; 
alpha_search_step = 0.05; % increase/decrease using stepsize of 5 percent


beta_LmaxNum = zeros(beta_num, 1);
beta_ROI_LmaxNum = zeros(beta_num, 1);
betaGrid = zeros(beta_num, 1);   % record of beta values searched


alpha_LmaxNum = zeros(alpha_num, 1);
alpha_ROI_LmaxNum = zeros(alpha_num, 1);
alphaGrid = zeros(alpha_num, 1); % record of alpha values searched


% Step 2.2 : Grid search along beta
np_temp = noiseParameter;

h = waitbar(0, 'Performing grid search wrt beta, Please wait...');
for m = 1 : beta_num
    np_temp(3) = noiseParameter(3) * (1 + (m - beta_num_half + 1) * beta_search_step);
    betaGrid(m) = np_temp(3); % save beta
    for n = 1 : optimImageStackLen
        fprintf('n = %f\n', n);
        [y, x] = alphaBeta(imGStack(:, :, n), imMaxStack(:, :, n), imMinStack(:, :, n), cStack(n).cands, np_temp, 0);
        temp_len = length(x)
        beta_LmaxNum(m) = beta_LmaxNum(m) + temp_len;
        for s = 1 : temp_len
            if (bw(round(y(s)), round(x(s))) == 1)
                beta_ROI_LmaxNum(m) = beta_ROI_LmaxNum(m) + 1;
            end
        end
    end
    waitbar(m / beta_num, h)
end
close(h);

Lmax_ratio = zeros(beta_num, 1)
for i = 1 : beta_num
    temp1 = beta_ROI_LmaxNum(i) / roi_area;
    temp2 = (beta_LmaxNum(i) - beta_ROI_LmaxNum(i)) / non_roi_area;
    Lmax_ratio(i) = temp1 / (temp1 + temp2);
end
%Lmax_ratio = Lmax_ratio * 

[temp_dif, temp_index] = min(abs(Lmax_ratio - confidenceP)); % Find the beta that gives the ratio closest to confidenceP

np_temp(3) = betaGrid(temp_index);
noiseParameter(3) = betaGrid(temp_index); % set the beta to the value found
% Step 2.3: Grid search along alpha

% surf(alphaBetaGrid);


% Step 3: --------------------------------------------------------
h = waitbar(0, 'Performing grid search wrt alpha, Please wait...');
for n = 1 : alpha_num
    np_temp(2) = noiseParameter(2) * (1 + (n - alpha_num_half + 1) * alpha_search_step);
    alphaGrid(n) = np_temp(2); % save beta
    for k = 1 : optimImageStackLen
        fprintf('k = %f\n', k);
        [y, x] = alphaBeta(imGStack(:, :, k), imMaxStack(:, :, k), imMinStack(:, :, k), cStack(k).cands, np_temp, 0);

        temp_len = length(x);
        alpha_LmaxNum(n) = alpha_LmaxNum(n) + temp_len;
        for s = 1 : temp_len
            if (bw(round(y(s)), round(x(s))) == 1)
                alpha_ROI_LmaxNum(n) = alpha_ROI_LmaxNum(n) + 1;
            end
        end
    end
    waitbar(n / alpha_num, h)
end
close(h);

%Lmax_ratio = alpha_ROI_LmaxNum ./ alpha_LmaxNum;
Lmax_ratio = zeros(alpha_num, 1)
for i = 1 : alpha_num
    temp1 = alpha_ROI_LmaxNum(i) / roi_area;
    temp2 = (alpha_LmaxNum(i) - alpha_ROI_LmaxNum(i)) / non_roi_area;
    Lmax_ratio(i) = temp1 / (temp1 + temp2);
end

[temp_dif, temp_index] = min(abs(Lmax_ratio - confidenceP)); % Find the beta that gives the ratio closest to confidenceP
noiseParameter(2) = alphaGrid(temp_index); % set the beta to the value found
actualP = Lmax_ratio(temp_index);

h = msgbox('Optimization of alpha & beta is done.');
pause(2);
close(h);

%--------------------------------------------------------------------------

function [y, x, cands]= alphaBeta(IG, Imax, Imin, cands, noiseParam, enhTriang, userROIbw)
% This is s simpifiled function for the statistical detection of speckles.
% It is written to avoid redundancy in speckle detection calculation. 

if nargin == 6
    userROIbw = []; % Set the optional value userROI to 0
end

if ~isempty(userROIbw)
    % Mask Imax
    Imax = Imax .* userROIbw;
end

% analyze speckles - validate, locmax, locmin...
[Imax, cands] = fsmPrepTestLocalMaxima(IG, Imax, cands, noiseParam, IG);  

% find the coordinates/positions of the local maxima after selecting only the significant local maxima/speckles
[y, x] = find(ne(Imax,0));

