L = zeros(201);
L(:,101) = 1;
Lg = imgaussfilt(L,2,'FilterSize',25);

box = zeros(201);
box(101+(-3:3),:) = 1;
Lg_box = imgaussfilt(L.*box,2,'FilterSize',25);

delta = zeros(201);
delta(101,101) = 1;
Dg = imgaussfilt(delta,3,'FilterSize',25);

testScales = 1/2/pi./chebpts(17,[1 4]);
F = orientationSpace.steerableScaleFilter(testScales,201);
R = F*(Lg+Lg.');
A = squeeze(R.getArraySpace);
B = squeeze(real(A(101,90,:,:)));
C = interpft(B,360);
cf = chebfun(C.',[1 4]);
[~,maxi] = max(cf);


%% Steerable scale filter

queryScales = 1:0.1:4;
queryOutput = zeros(size(queryScales));
parfor ii=1:length(queryScales)
    Lg = imgaussfilt(L,queryScales(ii),'FilterSize',25);
    R = F*Lg;
    A = squeeze(R.getArraySpace);
    B = squeeze(real(A(101,101,:,:)));
    cf = chebfun(B(1,:).',[1 4]);
    [~,queryOutput(ii)] = max(cf);
end


queryScales = 1:0.1:4;
queryOutput = zeros(size(queryScales));
parfor ii=1:length(queryScales)
    Lg_box = imgaussfilt(L.*box,queryScales(ii),'FilterSize',25);
    R = F*Lg_box;
    A = squeeze(R.getArraySpace);
    B = squeeze(real(A(101,101,:,:)));
    cf = chebfun(B(1,:).',[1 4]);
    [~,queryOutput(ii)] = max(cf);
end

queryScales = 1:0.1:4;
queryOutput = zeros(size(queryScales));
parfor ii=1:length(queryScales)
    Dg = imgaussfilt(delta,queryScales(ii),'FilterSize',25);
    R = F*Dg;
    A = squeeze(R.getArraySpace);
    B = squeeze(real(A(101,101,:,:)));
    cf = chebfun(B(1,:).',[1 4]);
    [~,queryOutput(ii)] = max(cf);
end

%% Ideal scale filter

sF = squeeze(orientationSpace.scaleFilter(chebpts(17,[1 4]),0));

queryScales = 1:0.1:4;
queryOutput = zeros(size(queryScales));
for ii=1:length(queryScales)
    Lg = imgaussfilt(L,queryScales(ii),'FilterSize',25);
    out = ifft2(bsxfun(@times,fft2(Lg),real(sF)));
    B = squeeze(real(out(101,101,:)));
    cf = chebfun(B,[1 4]);
    [~,queryOutput(ii)] = max(cf);
end


queryScales = 1:0.1:4;
queryOutput = zeros(size(queryScales));
for ii=1:length(queryScales)
    Dg = imgaussfilt(delta,queryScales(ii),'FilterSize',25);
    out = ifft2(bsxfun(@times,fft2(Dg),real(sF)));
    B = squeeze(real(out(101,101,:)));
    cf = chebfun(B,[1 4]);
    [~,queryOutput(ii)] = max(cf);
end