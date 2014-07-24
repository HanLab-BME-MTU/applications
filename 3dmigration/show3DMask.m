function show3DMask(maskIn)

fv = isosurface(maskIn);
patch(fv,'EdgeColor','none','FaceAlpha',.15)
[M,N,P] = size(maskIn);
xlim([0 N])
ylim([0 M])
zlim([0 P])
view(3);
axis equal




