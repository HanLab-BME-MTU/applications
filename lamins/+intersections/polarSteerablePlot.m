function polarSteerablePlot(steerable)
    % create image figure
    imFig = figure;
%     imshow(steerable.image,[]);
    imshowpair(steerable.image,steerable.nms > 0);
    h = impoint(gca,ceil(size(steerable.image,2)/2),ceil(size(steerable.image,1)/2));
    
    % setup accessory figures
    polarFig = figure;
    barFig = figure;
    fftFig = figure;

    
    % plot initial angle plots
    lin2 = @(x) [x(:) ;x(:); x(1)];
    plotPolar([ceil(size(steerable.image,2)/2),ceil(size(steerable.image,1)/2)]);
    plotBar([ceil(size(steerable.image,2)/2),ceil(size(steerable.image,1)/2)]);
    
    % add callbacks
    addNewPositionCallback(h,@plotPolar);
    addNewPositionCallback(h,@plotBar);
    addNewPositionCallback(h,@plotFFT);
    
    function plotPolar(x)
        figure(polarFig);
        theta = (0:5:360)/180*pi; 
        rho = steerable.a( round(x(2)) , round(x(1)) , : );
        rho = lin2(rho)';
        polar( 2*pi - theta, rho);
        hold on;
        polar( 2*pi- steerable.theta( round(x(2)) , round(x(1)) ) , steerable.res( round(x(2)) , round(x(1)) ),'ro');
        hold off;
%         disp(steerable.t( round(x(2)) , round(x(1)) ));
%         figure(imFig);
    end
    function plotBar(x)
        figure(barFig);
        theta = (0:5:360); 
        rho = steerable.a( round(x(2)) , round(x(1)) , : );
        rho = lin2(rho)';
        plot( 360 - theta, rho);
        set(gca,'XTick',0:45:360)
        hold on;
%         disp(steerable.t( round(x(2)) , round(x(1)) )/pi*180);
        plot( 360-mod(steerable.theta( round(x(2)) , round(x(1)) )/pi*180,360) , steerable.res( round(x(2)) , round(x(1)) ),'ro');
        hold off;
%         figure(imFig);
    end
    function plotFFT(x)
        figure(fftFig);
        theta = (0:5:360); 
        rho = steerable.a( round(x(2)) , round(x(1)) , : );
%         rho = [rho(:) ; rho(:) ]';
        rho = rho(:);
        fft_rho = fft(rho-mean(rho));
        mag = abs(fftshift(fft_rho));
        plot(mag,'b');
%         hold on;
%         phase = angle(fftshift(fft_rho))/pi;
%         plot(phase .* max(mag(:)),'r');
%         hold off;
%         set(gca,'XTick',0:45:360)
%         hold on;
%         disp(steerable.t( round(x(2)) , round(x(1)) )/pi*180);
%         plot( 360-mod(steerable.theta( round(x(2)) , round(x(1)) )/pi*180,360) , steerable.res( round(x(2)) , round(x(1)) ),'ro');
%         hold off;
%         figure(imFig);
    end
end