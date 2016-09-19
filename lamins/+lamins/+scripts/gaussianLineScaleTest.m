r = zeros(33);
F2 = OrientationSpaceRidgeFilter(1./2/pi./s,[],1,'none');
h = figure;
% for halflength = [0 10 25 50];
    for i=1:33;
        halflength = ceil(s(i)*3);
    %     G = repmat(normpdf(-49:50,0,s(i)),100,1);
    %     G = bsxfun(@times,normpdf(-49:50,0,s(i)),normpdf(-49:50,0,s(i)*10)');  
    %     G((-49:50)+50,:) = repmat(normpdf(-49:50,0,s(i)),100,1);
        G = zeros(201);
        G((-halflength:halflength)+101,50) = 1;
        G = imgaussfilt(G,s(i));
%         figure; imshow(G,[]);
        F2G = F2*G;
        F2Gas = F2G.getArraySpace();
        A = interpft(squeeze(F2Gas(101,50,:,1:10:33)),360);
        r(:,i) = squeeze(F2Gas(101,50,1,:))';
    end
    E = real(F2.getEnergy);
%     r = bsxfun(@times,r,1./real(E(:,1)).^1);
    cfs = chebfun(r,[1 8]);
    [maxr,maxs] = max(cfs);
    
    figure;
    plot(cfs)
    hold on;
    plot(maxs,maxr,'ko')
    xlabel('Scale (px)');
    ylabel('Response');
    title(sprintf('Object Length = %g',2*halflength+1));
    
    figure(h);
%     plot(s,maxs,'-')
    xlabel('Object Scale (px)')
    ylabel('Detected Scale (px)');
    maxcf = chebfun(maxs',[1 8]);
    hold on;
    plot(maxcf,'-o');
    hold on;
% end

figure(h)
plot([1 8],[1 8],'k--')
grid on
axis square

%% Show G fft
figure; imshow(abs(fftshift(fft2(G))),[]);