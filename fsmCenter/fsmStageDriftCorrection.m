function T = stageDriftCorrection(inputFileList, sigmaPSF, bitDepth, maxDrift, showResult)
% stageDriftCorrection returns an array containing all drifts between each
% consecutive pair of images.
%
% SYNOPSIS T = stageDriftCorrection(...)
%
% INPUT    inputFileList: cell array of all image filenames (including full
%                         path name).
%
%          sigmaPSF: hald-with of the point spread function, i.e. standard
%          deviation of a Gaussian model PSF (In pixel).
%
%          bitDepth: camera bit depth.
%
%          maxDrift: maximum drift allowed (in pixel).
%
%          showResult: show maximum projection 
%
% OUTPUT   T: an array of the same size of inputFileList minus 1 containing
%          all 2D drifts between images.
%
% Sylvain Berlemont, December 12th, 2009

T = [];

if nargin < 1 || isempty(inputFileList)
   [filename, pathname] = uigetfile({'*.tif';'*.jpg';'*.png';'*.*'}, ...
       'Select First Image');
   
   if ~ischar(filename) || ~ischar(pathname)
       return;
   end
   
   inputFileList = getFileStackNames([filename filesep pathname]);
else
    isValid = 1;
    for i = 1:numel(inputFileList)
        isValid = isValid && exist(inputFileList{i}, 'file');
    end
    if ~isValid
        error('Invalid input files.');
    end
end

if numel(inputFileList) < 2
    error('Number of images is less than 2.');
end

n = numel(inputFileList);
pts = cell(n, 1);
T = zeros(n - 1, 2);

%
% Step 1: Make calibration
%

[I0, sDN, GaussRatio] = fsmCalcNoiseParam(inputFileList{1}, bitDepth, sigmaPSF, []);
noiseParam = [1.96 / GaussRatio, sDN, 0, I0];

%
% Subpixel spot detection
%

h = waitbar(0, 'Subpixel spot detection');

for i=1:n
    % Load image
    I = double(imread(inputFileList{i}));
    
    % Normalized it
    In = I / (2^bitDepth-1);
    
    % Prepare the image for the analysis
    If = fsmPrepPrepareImage(In, 1, [1 1 0 0; 0 0 size(I)], sigmaPSF);
    
    % Statistically test the local maxima to extract (significant) speckles
    [~, cands] = fsmPrepMainSecondarySpeckles(If, 0, [], noiseParam, [1 0]);

    % Get speckle infos
    status = vertcat(cands(:).status);
    C = vertcat(cands(status == 1).Lmax);
    mu = mean(vertcat(cands(status == 1).IBkg)) * bitDepth;
    locMax = sub2ind(size(I), C(:, 1), C(:, 2));
    
    if ~numel(locMax)
        close(h);
        error('Frame %d does not contain any point.', i);
    end
    
    % Subpixel detection
    estimates = fitMixModel(I, [C(:, 1) C(:, 2)], sigmaPSF, I(locMax), mu);
    
    pts{i} = estimates(:,1:2);
    
    waitbar(i / n, h);
end

close(h);

%
% Step 3: Point registration
%

numIter = 10;
tol = 1e-4;

h = waitbar(0, 'Point registration');

for i = 1:n-1
    % Discard every point that doesn't contain a point in the next frame
    % within its vincinity,
    D = createDistanceMatrix(pts{i}, pts{i+1});
    D = D < maxDrift;
    X1 = pts{i}(sum(D, 2) ~= 0, 1:2);
    X2 = pts{i+1}(sum(D, 1) ~= 0, 1:2);
    n1 = size(X1, 1);
    n2 = size(X2, 1);
    
    if ~(n1 && n2)
        close(h);
        error('Unable to register frames %i-%i', i, i+1);
    end
    
    X1 = [X1 zeros(n1, 1)];
    X2 = [X2 zeros(n2, 1)];
    
    if n1 >= n2
        Ti = computeICP(X1, X2, numIter, tol);
    else
        Ti = -computeICP(X2, X1, numIter, tol);
    end
    
    T(i, :) = Ti(1:2);
    
    waitbar(i / (n-1), h);
end

close(h);

%
% Step 4: Register image for debug purposes
%

if showResult
    sumT = cumsum(T);
    maxX = ceil(max(abs(sumT(:, 2))));
    maxY = ceil(max(abs(sumT(:, 1))));
    I = double(imread(inputFileList{1}));
    I = padarray(I, [maxY, maxX]);
    %[path, ~, no] = getFilenameBody(inputFileList{1});
    %save([path filesep 'REGISTRED_' no '.txt'], 'I', '-ASCII', '-double');
    
    projI = I;
    projR = zeros(size(I));
    
    for i = 2:n
        % Read image.
        I = double(imread(inputFileList{i}));
        I = padarray(I, [maxY, maxX]);
        Tr = maketform('affine', [1 0 0; 0 1 0; fliplr(sumT(i-1, :)) 1]);
        R = imtransform(I, Tr, 'bicubic', 'XData', [1 size(I, 2)], ...
            'YData', [1 size(I, 1)]);
        %[path, ~, no] = getFilenameBody(inputFileList{i});
        %save([path filesep 'REGISTRED_' no '.txt'], 'R', '-ASCII', '-double');
        
        projI = projI + I;
        projR = projR + R;
    end
    
    colormap('jet');
    subplot(1, 2, 1); imagesc(projI); title('Original stack projection.');
    subplot(1, 2, 2); imagesc(projR); title('Registered stack projection.');
end