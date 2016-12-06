%% Parameters
pixelWidth = 2*pi/201;
% Standard deviation of gaussian in object space
obj_sigma = 2;
% Standard deviation of gaussian in frequency space
sigma_hat = 1/obj_sigma;
% Area under a gaussian from -Inf to Inf
gaussNorm = sqrt(2*pi)*sigma_hat;
% Half length
halfLength = 5;
% Sum of pixel values in object space
imgSum = halfLength*2+1;
format long;
theta = (-100:100)/201*2*pi;


%% Setup Line
L = zeros(201);
L((-halfLength:halfLength)+101,101) = 1;
Lg = imgaussfilt(L,2,'FilterSize',25);
L_hat = fftshift(fft2(ifftshift(L)));
Lg_hat = fftshift(fft2(ifftshift(Lg)));

%% Along the orientation axis
gaussSum = sum(real(Lg_hat(101,:)));
% sum is the theoretical gaussian sum divided by pixel width by object sum
gaussSumExpression = sum(exp(-theta.^2/2/sigma_hat.^2)*imgSum);
disp('Along orientation axis:');
disp(gaussSum);
disp(gaussNorm/pixelWidth*imgSum);
disp(sum(gaussSumExpression));
assert(abs(gaussSum - gaussNorm/pixelWidth*imgSum) < 1e-6);


%% Perpendicular to the orientation axis
diricSum = sum(real(Lg_hat(:,101)));
% sqrt(2*pi)/sigma_hat is the area of the inverse gaussian
% gaussNorm/pixelWidth is the integral of the original gaussian divided by
% the image sum since it's absorbed by diric
diricSumExpression = sum(exp(-theta.^2/2/sigma_hat.^2)*imgSum.*diric(theta,imgSum));
k = -halfLength:halfLength;
diricSymbolic = sum(exp(-k.^2*sigma_hat.^2/2))*gaussNorm/pixelWidth;
disp('Perpendicular to orientation axis:');
disp(diricSum);
disp(diricSymbolic);
% disp(sqrt(2*pi)/sigma_hat*gaussNorm/pixelWidth);
disp(diricSumExpression);
assert(abs(diricSum - diricSymbolic) < 1e-7);

%% Reset format
format
