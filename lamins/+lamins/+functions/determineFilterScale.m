function [scale] = determineFilterScale(template,start,finish,TOL)

    if(nargin < 2)
        start = 1;
    end
    if(nargin < 3)
        finish = 8;
    end
    if(nargin < 4)
        TOL = 1e-6;
    end
    if(nargin < 1 || isempty(template))
        disp('Generating template');
        templateScale = rand(1)*(finish-start)+start;
        fprintf('template scale = %g\n',templateScale);
        template = zeros(100);
        template(50,50) = 1;
        template = orientationSpace.kernel(1/2/pi/templateScale,[],8,0,100).*fft2(template);
        template = ifft2(template);
    end
    center.r = floor(size(template,1)/2);
    center.c = floor(size(template,2)/2);
    npts = 513;
    s = chebpts(npts,[start finish]);
    iters = log2(length(s)-1)-1;
    F2 = OrientationSpaceRidgeFilter(1./2/pi./s,[],8,'energy');
    F2([1 npts]).setupFilter(size(template));
    F2R([1 npts],1) = F2([1 npts])*template;
    F2([1 npts]).clearTransients;
    scale = zeros(1,iters+1);
    figure;
    for i=1:iters+1
    % for p=2.^(iters:-1:0)
        p = 2.^(iters+1-i);
        idx = (1+p):p*2:npts;
        F2(idx).setupFilter(size(template));
        F2R(idx) = F2(idx)*template;
        F2(idx).clearTransients;
        idx = 1:p:npts;

        F2Ras = F2R(idx).getArraySpace;
        r = squeeze(F2Ras(center.r,center.c,1,:));
        % r = r./real(F2(idx).getEnergy).^0.6105;
        cf = chebfun(r,[start finish]);
        h = plot(cf,'-o','MarkerSize',i*3+3);
        hold on
        [~,scale(i)] = max(cf)
        plot([scale(i) scale(i)],ylim,'Color',h.Color);
        pause(0.1);
        if(i > 1 && abs(scale(i) - scale(i-1)) < TOL)
            break;
        end
        pause(1);
    end
    scale = scale(1:i);

    clear F2R

    figure; plot(scale);

    if(nargin < 1)
        fprintf('Scale Error %0.2g%%\n',(1-templateScale/scale(end))*100);
    end

end


    % F = OrientationSpaceRidgeFilter(1./2/pi./3.1428,[],8); F.setupFilter([ 100 100]); F.F = F.F*3.1428;
    % FR = F*template;
    % FR.a(50,50,1)
