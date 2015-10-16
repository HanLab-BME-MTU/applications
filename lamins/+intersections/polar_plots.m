%% + Cross

plusBW.image = zeros(1024);
plusBW.image(:,512) = 1;
plusBW.image(512,:) = 1;
plusBW.image = imfilter(plusBW.image,fspecial('gaussian',25,5));
% plusBW.image = imnoise(plusBW.image,'gaussian');
[plusBW.res,plusBW.theta,plusBW.nms,plusBW.a] = steerableDetector(plusBW.image,4,5);

%% X cross
xBW.image = eye(1024);
xBW.image = double(xBW.image | fliplr(xBW.image));
xBW.image = imfilter(xBW.image,fspecial('gaussian',25,5));
xBW.image = imnoise(xBW.image,'gaussian',0,1e-4);
[xBW.res,xBW.theta,xBW.nms,xBW.a] = steerableDetector(xBW.image,4,5);

%% > cross
gtBW.image = zeros(1024);
gtBW.image(:,512) = 1;
gtBW.image = double(gtBW.image | eye(1024));
gtBW.image = imfilter(gtBW.image,fspecial('gaussian',25,5));
[gtBW.res,gtBW.theta,gtBW.nms,gtBW.a] = steerableDetector(gtBW.image,4,5);

%% actual data
L.image = double(images{2}(3,15).image);
[L.res,L.theta,L.nms,L.a] = steerableDetector(L.image,4,5);

%% Show images
figure;
subplot(1,3,1);
imshow(plusBW.image,[]);
subplot(1,3,2);
imshow(xBW.image,[]);
subplot(1,3,3);
imshow(gtBW.image,[]);

%% + Cross test
polarSteerablePlot(plusBW);

%% X cross test
polarSteerablePlot(xBW);

%% > cross test
polarSteerablePlot(gtBW);

%% actual data test
polarSteerablePlot(L);