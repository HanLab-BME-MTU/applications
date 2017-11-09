function [ output_args ] = interactiveVanGinkel( I, gamma )
%interactiveVanGinkel Interactive vanGinkel

if(nargin < 2)
    gamma = 0.4;
end

If = fftshift(fft2(I));

figure;
him = imshow(mat2gray(abs(If)).^gamma);
h = lamins.simsim.imRadialLine(gca,[512 512; 512 562]);
h2 = lamins.simsim.imThetaOffsetLine(h,gca,[512 512; 562 512]);

dim.x = xlim;
dim.y = ylim;

center = [mean(dim.x) mean(dim.y)];
hf = [];

while(isvalid(h))
    wait(h);
    if(isa(h,'imline') && isvalid(h))
        pos{1} = getPosition(h);
    else
        break;
    end
    [theta(1,1),rho(1,1)] = cart2pol(pos{1}(1,1)-center(1),pos{1}(1,2)-center(2));
    [theta(1,2),rho(1,2)] = cart2pol(pos{1}(2,1)-center(1),pos{1}(2,2)-center(2));
    
    pos{2} = getPosition(h2);
    [theta(2,1),rho(2,1)] = cart2pol(pos{2}(1,1)-center(1),pos{2}(1,2)-center(2));
    [theta(2,2),rho(2,2)] = cart2pol(pos{2}(2,1)-center(1),pos{2}(2,2)-center(2));
    
    N = size(I,1);
    
    f_c = mean(max(rho,[],2))/N;
    b_f = mean(abs(diff(rho,1,2)))/N;
    K = 2*pi/mean(abs(diff(theta)));
    K = ceil(max(K,1/(1-1/K)));
    angle = mean(theta(:));
    
    F = steerableVanGinkelKernel(f_c,b_f,K,angle,N/2);
    
%     figure;
    out = repmat(mat2gray(abs(If)).^gamma,[1 1 3]);
%     out(:,:,2) = abs(imag(F)).*mat2gray(imag(F)).*abs(If);
    out(:,:,3) = mat2gray(real(F).*abs(If)).^gamma;
%     out = mat2gray(out).^gamma;
    set(him,'CData',out);
    
    fprintf('f_c: %f\n',f_c);
    fprintf('b_f: %f\n',b_f);
    fprintf('K: %f\n',K);
    fprintf('angle: %f\n',angle);
%     imshow(out,[]);
    
    if(isempty(hf) || ~ishandle(hf) || ~isvalid(hf))
        hf = figure;
    else
        hf = figure(hf);
    end
    I_filtered = apply_freq_filter(fftshift(If),fftshift(F));
    imshowpair(real(I_filtered),abs(imag(I_filtered)));
end

function R = apply_freq_filter(If,f)
    R = ifft2(If.*real(f)) + 1j*ifft2(If.*imag(f).* -1j);
end

end

