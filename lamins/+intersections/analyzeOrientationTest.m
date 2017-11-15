figure;
for i=0:0.1:12
Ig = intersections.drawGaussianLine([pi/12 -pi/12]*i,2,[51 51],30);
Ign = imnoise(Ig,'gaussian',0,1e-2);
R = F*Ign;
[mm,~,mv] = interpft_extrema(real(R.a(51,51,:)),3);
mm(mv < mean(real(R.a(51,51,:)),3)) = NaN;
plot(i/12*180,mm/2/pi*180,'.'); hold on;
end