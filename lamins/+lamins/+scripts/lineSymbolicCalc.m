%% Parameters
pixelWidth = 2*pi/201;
% Standard deviation of gaussian in object space
obj_sigma = 2;
% Standard deviation of gaussian in frequency space
sigma_hat = 1/obj_sigma;
% Area under a gaussian from -Inf to Inf
gaussNorm = sqrt(2*pi)*sigma_hat;
% Sum of pixel values in object space
imgSum = 31;


%% Setup Line
L = zeros(201);
L((-15:15)+101,101) = 1;
Lg = imgaussfilt(L,2,'FilterSize',25);
L_hat = fftshift(fft2(ifftshift(L)));
Lg_hat = fftshift(fft2(ifftshift(Lg)));

%% Along the orientation axis
gaussSum = sum(real(Lg_hat(101,:)));
% sum is the theoretical gaussian sum divided by pixel width by object sum
assert(abs(gaussSum - gaussNorm/pixelWidth*imgSum) < 1e-6);
sum(exp(-theta.^2/2/sigma.^2)*31);

%% Perpendicular to the orientation axis
diricSum = sum(real(Lg_hat(:,101)));
% sqrt(2*pi)/sigma_hat is the area of the inverse gaussian
% gaussNorm/pixelWidth is the integral of the original gaussian divided by
% the image sum since it's absorbed by diric
assert(abs(diricSum - sqrt(2*pi)/sigma_hat*gaussNorm/pixelWidth) < 1e-14);
sum(exp(-theta.^2/2/sigma.^2)*31.*diric(theta,31));

