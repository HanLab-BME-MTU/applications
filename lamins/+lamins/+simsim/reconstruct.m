function [ I6 ] = reconstruct( I3,PSF )
%reconstruct Perform SIM reconstruction

    OTF = fftshift(fft2(PSF,1024,1024));
    MTF = abs(OTF);

    % Apply FFT to structured illuminated images
    I3f = cellfun(@(x) fftshift(fft2(x)),I3,'Unif',false);
    % Stack the images
    % I3fs = cat(3,I3f{:});
    I3fs = cell2mat(shiftdim(I3f,-2));

    % Rotate the matrix dimensions around to 3 x n , where N = X x Y
    I4 = shiftdim(I3fs,2);
    I4 = I4(:,:);

    % Construct systems of equations coefficients
    C = zeros(3);
    C(:,1) = 1;
    C(:,2) = 0.5*exp(1j*2*pi/3*[1 2 3]);
    C(:,3) = 0.5*exp(-1j*2*pi/3*[1 2 3]);

    % Solve systems of equations and rotate dimensions back
    I5 = shiftdim(reshape(C\I4,3,3,1024,1024),2);

    % Visualize
    % p.x = xlim; p.y = ylim;
    p.x = [-100 100]+513;
    p.y = [-100 100]+513;
    figure;
    for i=1:3;
    subplot(1,3,i);
    imshow(mat2gray(abs(I5(:,:,:,i))),[]); xlim(p.x); ylim(p.y);
    end


    % Reassign to higher frequencies
    I6 = I5;
    % freqshift = 1024/winding;
    N = repmat(OTF,[1 1 3 3]);
    doshift = @(I,angle_j,winding) fftshift(fft2(ifft2(fftshift(I)).*exp(1j.*Xr{angle_j}.*2*pi/winding)));
    for angle_j = 1:3
        I6(:,:,1,angle_j) = I6(:,:,1,angle_j).*MTF;
        I6(:,:,2,angle_j) = doshift(I5(:,:,2,angle_j).*MTF,angle_j,-winding);
        I6(:,:,3,angle_j) = doshift(I5(:,:,3,angle_j).*MTF,angle_j,winding);
        % I6(:,:,2) = circshift(I5(:,:,2).*OTF,-freqshift,2);
        % I6(:,:,3) = circshift(I5(:,:,3).*OTF,freqshift,2);


        N(:,:,2,angle_j) = doshift(N(:,:,2,angle_j),angle_j,-winding);
        N(:,:,3,angle_j) = doshift(N(:,:,3,angle_j),angle_j,winding);
    end

    N = abs(N);
    N = sum(sum(N.^2,4),3)+1e-1;

    I6 = I6./repmat(N,[1 1 3 3]);


end

