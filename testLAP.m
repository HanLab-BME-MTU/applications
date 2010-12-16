function testLAP(movieData,thE,thA,thP)

% thE: [0, +inf)
% thA: [0, pi]
% thP: [0, 1]

imagePath = fullfile(movieData.imageDirectory, movieData.channelDirectory{1});
imageFiles = dir([imagePath filesep '*.tif']);
ima = imread(fullfile(imagePath, imageFiles(1).name));

load(fullfile(movieData.particleDetection.directory, ...
    movieData.particleDetection.filename));

X = [featuresInfo(1).xCoord, featuresInfo(1).yCoord];
ind = sub2ind(size(ima),X(:,2),X(:,1));
N = size(ind,1);

[~, T] = steerableFiltering(double(ima),2,2);

Y = [cos(T(ind)), sin(T(ind))];

pair = pcombs(1:N,false);

u0 = X(pair(:,1),:);
u1 = u0 + Y(pair(:,1),:);
v0 = X(pair(:,2),:);
v1 = v0 + Y(pair(:,2),:);

isValid = true(size(pair,1),1);

% euclidian distance [0...+inf]
dE = sqrt(sum((u0 - v0).^2,2));
isValid = isValid & dE <= thE;

% angle between u and v
% dot = abs(sum((u1 - u0) .* (v1 - v0),2));
% dot(dot > 1) = 1;
% dot(dot < -1) = -1;
% dA = acos(dot);
% isValid = isValid & dA <= thA;
% 
% mean distance of u1 and v1 projected on the line (u0,v0)
dP1 = sqrt(1 - (sum((u1-u0) .* (v0-u0),2) ./ dE).^2);
dP2 = sqrt(1 - (sum((v1-v0) .* (v0-u0),2) ./ dE).^2);
dP = dP1 .* dP2;
isValid = isValid & dP <= thP;
% 
% cost = exp(- (dE .* (1/pi) .* dA .* dP));

cost = 1./ dE;

% Build cost matrix
i = pair(isValid,1);
j = pair(isValid,2);
c = cost(isValid);

% Populate the lower triangular part only
D = sparse(j, i, c, N, N, numel(c));

% Compute Maximum Weight Matching
M = maxWeightMatching(D);


% Display result
imshow(ima,[]);
hold on;

B = zeros(N);
ind = sub2ind([N N],M(:,1),M(:,2));
B(ind) = 1;

line(X(:,1),X(:,2),'LineStyle','none', 'Marker', '.', 'Color', 'g');
gplot(B,X,'r');
quiver(X(:,1),X(:,2),Y(:,1),Y(:,2),0,'b');

end