deltaAngle = 360/197;
endAngle = 360;
angles = 0:deltaAngle:endAngle;
I = Lamins(2).cellReader{3,1,15};
mask = Lamins(2).getNuclearMask(3,1,15);
I = im2double(imadjust(I,stretchlim(I,0)));
bins = -4:0.1:4;

orientations = cell(length(angles),1);
maxima = zeros(1,length(angles));

parfor ii=1:length(angles)
    J = imrotate(I,angles(ii));
    Jmask = imrotate(mask,angles(ii));
    [res, orientations{ii}, nms] = steerableDetector(J,4,5);
    orientations{ii} = orientations{ii}.*sign(nms);
    orientations{ii} = orientations{ii}(Jmask & nms > 0);
    [~,maxima(ii)] = max(histc(orientations{ii}/pi*8,bins))
end

figure;
[res, theta, nms] = steerableDetector(I,4,5);
imshow((theta*180/pi+90).*(nms > 0).*mask-90,[]);
cm = [ [0 0 0]; hsv];
colormap(cm);
thetaC = (theta*180/pi+90).*(nms > 0).*mask;
thetaC = ceil(thetaC);
% imwrite(thetaC,[ [0 0 0]; hsv(180) ],'Fig_04_Unrotated_theta.png')

figure;
histogram(theta/pi*180);

figure;
plot(angles,bins(maxima)/8*180,'-g');
hold on;
plot(angles,bins(maxima)/8*180,'o');
hold off;
