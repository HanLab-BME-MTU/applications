
%% Construct test image
I = zeros(1024);

% Draw two lines
out = bresenham([1 450],[1024 550]);
I(sub2ind([1024 1024],out(:,1),out(:,2))) = 1;
% out = bresenham([1 900],[1024 124]);
out = bresenham([450 1],[550 1024]);
I(sub2ind([1024 1024],out(:,1),out(:,2))) = 1;

% Delete 80% of the points from the line
linepoints = find(I);
linepoints = linepoints(randperm(length(linepoints),round(length(linepoints)*80/100)));
I(linepoints) = 0;

% Also scatter 3000 random points throughout
I(randi(numel(I),3000,1)) = 1;
figure;
imshow(I);


%% Setup PSF
psf_p.NA = 1.2;
% 31.6 nm
psf_p.pixel_size = sqrt(1000)/1000;
psf_p.lambda = 0.52;
PSF = psf2D(psf_p.pixel_size,psf_p.NA,psf_p.lambda);
% PSF = fspecial('gaussian',50,8);
OTF = fftshift(fft2(PSF,1024,1024));
MTF = abs(OTF);

% Apply a Gaussian 'PSF'
I2 = imfilter(I,PSF);

% Do structured illumination
winding = 16;
[X,Y] = meshgrid(1:1024);
theta = [1 2 3]*2*pi/3;
% X = sin(theta).*X - cos(theta).*Y;

I3 = cell(3,3);
L = cell(3,3);
Xr = cell(1,3);
Yr = cell(1,3);
for angle_j=1:3
    Xr{angle_j} = cos(theta(angle_j)).*X - sin(theta(angle_j)).*Y;
    Yr{angle_j} = sin(theta(angle_j)).*X + cos(theta(angle_j)).*Y;
    for phase_i=1:3
        L{phase_i,angle_j} = (1+cos(Xr{angle_j}*2*pi/winding+phase_i*2*pi/3))/2;
        I3{phase_i,angle_j} = imfilter(I.*L{phase_i,angle_j},PSF);
    end;
end

figure;
imshowpair(abs(fftshift(fft2(L{1}))),MTF);

% Apply FFT to structured illuminated images
I3f = cellfun(@(x) fftshift(fft2(x)),I3,'Unif',false);
% Stack the images
% I3fs = cat(3,I3f{:});
I3fs = cell2mat(shiftdim(I3f,-2));

% Rotate the matrix dimensions around to 3 x n , where N = X x Y
I4 = shiftdim(I3fs,2);
I4 = I4(:,:);

% Construct systems of equations coefficients
C = zeros(3);
C(:,1) = 1;
C(:,2) = 0.5*exp(1j*2*pi/3*[1 2 3]);
C(:,3) = 0.5*exp(-1j*2*pi/3*[1 2 3]);

% Solve systems of equations and rotate dimensions back
I5 = shiftdim(reshape(C\I4,3,3,1024,1024),2);

% Visualize
% p.x = xlim; p.y = ylim;
p.x = [-100 100]+513;
p.y = [-100 100]+513;
figure;
for i=1:3;
subplot(1,3,i);
imshow(mat2gray(abs(I5(:,:,:,i))),[]); xlim(p.x); ylim(p.y);
end


% Reassign to higher frequencies
I6 = I5;
% freqshift = 1024/winding;
N = repmat(OTF,[1 1 3 3]);
doshift = @(I,angle_j,winding) fftshift(fft2(ifft2(fftshift(I)).*exp(1j.*Xr{angle_j}.*2*pi/winding)));
for angle_j = 1:3
    I6(:,:,1,angle_j) = I6(:,:,1,angle_j).*MTF;
    I6(:,:,2,angle_j) = doshift(I5(:,:,2,angle_j).*MTF,angle_j,-winding);
    I6(:,:,3,angle_j) = doshift(I5(:,:,3,angle_j).*MTF,angle_j,winding);
    % I6(:,:,2) = circshift(I5(:,:,2).*OTF,-freqshift,2);
    % I6(:,:,3) = circshift(I5(:,:,3).*OTF,freqshift,2);

    
    N(:,:,2,angle_j) = doshift(N(:,:,2,angle_j),angle_j,-winding);
    N(:,:,3,angle_j) = doshift(N(:,:,3,angle_j),angle_j,winding);
end

N = abs(N);
N = sum(sum(N.^2,4),3)+1e-1;

I6 = I6./repmat(N,[1 1 3 3]);


%% Visualize reconstruction in frequency space
figure;
imshow(mat2gray(abs(sum(I6,4))).^0.4,[]); xlim(p.x); ylim(p.y);
figure;
imshow(mat2gray(abs(squeeze(sum(I6,3)))).^0.4); xlim(p.x); ylim(p.y);

%% Visualize reconstruction in real space
I8 = ifft2(fftshift(sum(sum(I6,4),3)));
figure;
imshow(I8,[])
figure;
imshowpair(I2,abs(I8));
