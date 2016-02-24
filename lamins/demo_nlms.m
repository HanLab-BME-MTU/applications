function out = demo_nlms
[s,rho,I] = intersections.synthesizeIntersection(0:pi/6:pi-pi/6);
[v.res, v.theta, v.nms, v.a] = steerableOrientationSpaceFilter(I,0.05,0.04,6,256);
% [v.res, v.theta, v.nms, v.a] = steerableDetector(I,4,5,256);
nlms = nonLocalMaximaSuppression(real(v.a));
T = thresholdOtsu(nlms(nlms > 0));
spy3d(nlms > T)
hold on;
surf(-ones(100),I,'EdgeColor','none');
zlim([-1 256]);
colormap gray
out = nlms;
end