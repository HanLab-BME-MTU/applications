s = chebpts(33,[1 4]);

%% 2D Gaussian Spot

r = zeros(33);
F2 = OrientationSpaceRidgeFilter(1./2/pi./s,[],0,'none');
for i=1:33;
    G = imgaussfilt(delta,s(i),'FilterSize',41);
    F2G = F2*G;
    F2Gas = F2G.getArraySpace();
    r(:,i) = squeeze(F2Gas(50,50,1,:))';
end

imagesc(r)
imagesc(bsxfun(@rdivide,r,max(r)))
imagesc(bsxfun(@rdivide,r,max(r)) == 1)
rsqrts = bsxfun(@times,r,1./real(F2.getEnergy).^(0.6105));
imagesc(bsxfun(@rdivide,rsqrts,max(rsqrts)))
imagesc(bsxfun(@rdivide,rsqrts,max(rsqrts)) == 1)

%% Determine Gaussian 2D Scale Demo

lamins.functions.determineGaussian2DScale

%% Gaussian Line

r = zeros(33);
F2 = OrientationSpaceRidgeFilter(1./2/pi./s,[],0,'none');
for i=1:33;
    G = repmat(normpdf(-49:50,0,s(i)),100,1);
%     G = bsxfun(@times,normpdf(-49:50,0,s(i)),normpdf(-49:50,0,s(i)*3)');  
%     G((-5:5)+50,:) = repmat(normpdf(-49:50,0,s(i)),11,1);
    F2G = F2*G;
    F2Gas = F2G.getArraySpace();
    r(:,i) = squeeze(F2Gas(50,50,1,:))';
end

imagesc(r)
imagesc(bsxfun(@rdivide,r,max(r)))
imagesc(bsxfun(@rdivide,r,max(r)) == 1)
% rsqrts = bsxfun(@times,r,1./real(F2.getEnergy).^(0.6105));
% imagesc(bsxfun(@rdivide,rsqrts,max(rsqrts)))
% imagesc(bsxfun(@rdivide,rsqrts,max(rsqrts)) == 1)

% lamins.functions.determineGaussianLineScale

%% Filter pattern

r = zeros(33);
F2 = OrientationSpaceRidgeFilter(1./2/pi./s,[],5,'none');
for i=1:33;
    G = zeros(100);
    G(50,50) = 1;
    G = orientationSpace.kernel(1/2/pi/s(i),[],8,0,100).*fft2(G);
    G = ifft2(G);
    F2G = F2*G;
    F2Gas = F2G.getArraySpace();
    r(:,i) = squeeze(F2Gas(50,50,1,:))';
end

imagesc(r)
imagesc(bsxfun(@rdivide,r,max(r)))
imagesc(bsxfun(@rdivide,r,max(r)) == 1)
E = real(F2.getEnergy);
rsqrts = bsxfun(@times,r,1./real(E(:,1)));
imagesc(bsxfun(@rdivide,rsqrts,max(rsqrts)))
imagesc(bsxfun(@rdivide,rsqrts,max(rsqrts)) == 1)

%% 