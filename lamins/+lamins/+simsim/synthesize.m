function [ I2, PSF ] = synthesize( )
%synthesize Create a synthetic image

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
psf_p.lambda = 0.561;
PSF = psf2D(psf_p.pixel_size,psf_p.NA,psf_p.lambda);
% PSF = fspecial('gaussian',50,8);
% OTF = fftshift(fft2(PSF,1024,1024));
% MTF = abs(OTF);

% Apply a Gaussian 'PSF'
I2 = imfilter(I,PSF);

end

