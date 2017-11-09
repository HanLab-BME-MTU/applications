I = intersections.drawTwoLines([0 pi/3]);
figure; imshow(I,[]);
Ig = imgaussfilt(I,2,'FilterSize',25);
Ign = imnoise(Ig,'gaussian',0,1e-3)
figure; imshow(I,[]);
figure; imshow(Ign,[]);
R = F*Ign;
r = 51; c = 51;
K = 8:-0.1:1;
rho = R.getResponseAtOrderFTatPoint(r,c,K);
rhoh = fft(rho);
rhoh(1,:) = 0;
rhoh(abs(rhoh) < eps*1e3) = 0;
rho = ifft(rhoh);
out = interpft_extrema(rhoh,1,[],[],false);
% out = interpft_extrema(rho);
[out,out_events] = orientationSpace.diffusion.alignExtrema(out);
out = orientationSpace.diffusion.unwrapExtrema(out,out_events);

figure; plot(K,out/2/pi*180);