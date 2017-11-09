angles_degrees = 0:5:180;
angles = angles_degrees*pi/180;
x = [50 50];
lin2 = @(x) [x(:) ;x(:); x(1)];

for a = 1:length(angles);
    lineAngle = angles(a);
    basename = ['angle ' num2str(angles_degrees(a))];
    I = drawTwoLines(lineAngle);
    I = imfilter(I,fspecial('gaussian',10,2));
    I = imnoise(I,'gaussian',0,1e-4);
    imshow(I);
    f = getframe(gcf);
    imwrite(f.cdata,[basename '.png']);
    %pause;
    
    [s.res, s.theta, s.nms, s.a] = steerableDetector(double(I),4,2);
    imshow(s.nms,[]);
    f = getframe(gcf);
    imwrite(f.cdata,[basename '_nms.png']);
    %pause;
    
    theta = (0:5:360)/180*pi; 
    rho = s.a( round(x(2)) , round(x(1)) , : );
    rho = lin2(rho)';
    polar( 2*pi - theta, rho);
    hold on;
    polar( 2*pi- s.theta( round(x(2)) , round(x(1)) ) , s.res( round(x(2)) , round(x(1)) ),'ro');
    hold off;
    f = getframe(gcf);
    imwrite(f.cdata,[basename '_polar.png']);
    %pause;
   
   plot(theta,rho);
   xlabel('Orientation (radians)');
    f = getframe(gcf);
    imwrite(f.cdata,[basename '_plot.png']);
    %pause;
    
    reciprical = (-36:35)/2/pi/72;
    rho = rho(1:end-1);
    fft_rho = fft(rho-mean(rho));
    mag = abs(fftshift(fft_rho));
    plot(reciprical,mag,'b');
    xlabel('1/Orientation (1/radians)');
    ylabel('Magnitude');
    f = getframe(gcf);
    imwrite(f.cdata,[basename '_fft.png']);
    %pause;
    
    fft_angle = angle(fftshift(fft_rho));
    plot(reciprical,fft_angle,'r');
    xlabel('1/Orientation (1/radians)');
    ylabel('Phase angle (radians)');
    f = getframe(gcf);
    imwrite(f.cdata,[basename '_phase.png']);
    %pause;
end