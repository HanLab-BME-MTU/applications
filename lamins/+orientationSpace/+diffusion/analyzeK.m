if(~exist('I','var'))
    [~,hostname] = system('hostname');
    hostname = strtrim(hostname);

    switch(hostname)
        case 'FSM2UA220003Q'
            % HP Z820 Goldman workstation
            cd 'Z:\Takeshi\N-SIM\040715';
            MD = MovieData.load('MEFLB1-LACLB12-006_Reconstructed.nd2');
            I = MD.channels_(1).loadImage(1,11);
        case 'mkitti-jaqaman'
            % Laptop T440s
            load('C:\Users\Mark Kittisopikul\Documents\Data\Lamins\MEFLB1-LACLB12-006_Reconstructed_study\MEFLB1-LACLB12-006_Reconstructed\MEFLB1-LACLB12-006_Reconstructed.mat');
            MD.sanityCheck;
            I = MD.channels_(1).loadImage(1,10);
        otherwise
            % BioHPC
            cd ~/shortcuts/MEFLB1-LACLB12-006_Reconstructed/
            MD = MovieData.load('MEFLB1-LACLB12-006_Reconstructed.mat');
            I = MD.channels_(1).loadImage(1,10);
    end
end
% I = imread('example.tif');
F = OrientationSpaceFilter.constructByRadialOrder(1/2/pi./2,1,8,'none');
R = F*I;

% Obtain derivative responses
% a = ifft(bsxfun(@times,fft(shiftdim(real(R.a),2)),shiftdim(ifftshift(-8:8).*1i,1)));
% Rd = OrientationSpaceResponse(R.filter,shiftdim(a,1));
% a = ifft(bsxfun(@times,fft(shiftdim(real(R.a),2)),shiftdim(ifftshift(-8:8).*1i,1).^2));
% Rdd = OrientationSpaceResponse(R.filter,shiftdim(a,1));
Rd  =  R.getDerivativeResponse(1);
Rdd = Rd.getDerivativeResponse(1);

% Get Maxima
maxima = cell(1,8);
for i=1:8; maxima{i} = R.getResponseAtOrderFT(i,2).getRidgeOrientationLocalMaxima; end;

% Get Maxima, Minima of Derivative
maxima_d = cell(1,8);
minima_d = cell(1,8);
for i=1:8; [maxima_d{i},minima_d{i}] = Rd.getResponseAtOrderFT(i,2).getRidgeOrientationLocalMaxima; end;


% Make them all the same size
% for i=1:7;
%     maxima{i}(:,:,end+1:7) = NaN;
%     maxima_d{i}(:,:,end+1:8) = NaN;
%     minima_d{i}(:,:,end+1:8) = NaN;
% end;
% 
% % Row x Column x Extrema x K
% maxima = cat(4,maxima{:});
% maxima_d = cat(4,maxima_d{:});
% minima_d = cat(4,minima_d{:});

%% Construct interpolation matrix
tau = augknt([0 1],3);
knots = aptknt(tau,4);
colmat = spcol(knots,4,augknt([0 1],3),'sl');

%% Extract specific pixel info

r = 629; c = 373;

ms = max(cellfun('size',maxima,3));
m = NaN(ms,length(maxima));
for i=1:8
    m(1:size(maxima{i}(r,c,:),3),i) = squeeze(maxima{i}(r,c,:));
end
m = orientationSpace.diffusion.alignExtrema(m);
ms = max(cellfun('size',maxima_d,3));
maxd = NaN(ms,length(maxima_d));
for i=1:8
    maxd(1:size(maxima_d{i}(r,c,:),3),i) = squeeze(maxima_d{i}(r,c,:));
end
maxd = orientationSpace.diffusion.alignExtrema(maxd);
ms = max(cellfun('size',minima_d,3));
mind = NaN(ms,length(minima_d));
for i=1:8
    mind(1:size(minima_d{i}(r,c,:),3),i) = squeeze(minima_d{i}(r,c,:));
end
mind = orientationSpace.diffusion.alignExtrema(mind);

%% Combine maxima into a single matrix

ms = max(cellfun('size',maxima,3));

for i=1:length(maxima)
    maxima{i}(:,:,end+1:ms) = NaN;
end
% R x C x LM x K
maxima = cat(4,maxima{:});
% LM x K x R x C
maxima = shiftdim(maxima,2);

parfor r=1:1024
    for cc=1:1024
        maxima(:,:,r,cc) = orientationSpace.diffusion.alignExtrema(maxima(:,:,r,cc),pi,true,false);
    end;
    fprintf('r = %d, c = %d\n',r,cc);
end;

%% Do rough maxima count

maxima_nan = isnan(maxima);
maxima_count = ms-sum(maxima_nan);
maxima_count = squeeze(shiftdim(maxima_count,2));
maxima_count_diff = diff(maxima_count,1,3) > 0;
loss_labeled = bsxfun(@times,maxima_count_diff,shiftdim(1:ms,-1));
figure;
for i=ms:-1:1
    loss_map= max(loss_labeled(:,:,1:i),[],3);
    imshow(loss_map,[0 ms]);
    colormap(gca,parula);
%     caxis([0 ms]);
    pause;
end

%% Find maxima that end

maxima_stop_pts = maxima_nan(:,1:7,:,:) & ~maxima_nan(:,2:8,:,:);
sz = size(maxima_stop_pts);
maxima_stop_pts = [false(sz(1),1,sz(3),sz(4)) maxima_stop_pts];

%% Find precise maxima stop point

zeroCross = cell(1024,1024,8);
parfor K=2:8
    derivs = orientationSpace.diffusion.orientationMaximaDerivatives(shiftdim(real(R.a),2),K,3,squeeze(maxima(:,K,:,:))*2);
    derivs = permute(derivs,[1 4 2 3]);
    for r=1:1024
        for c=1:1024
            s = maxima_stop_pts(:,K,r,c);
            p = fliplr(bsxfun(@rdivide,[maxima(s,K,r,c)*2 derivs(s,:,r,c)],factorial(0:3)));
            rho_d = Rd.getResponseAtOrderFTatPoint(r,c,K-1:0.01:K);
            rho_dd = Rdd.getResponseAtOrderFTatPoint(r,c,K-1:0.01:K);
            % Relative K
            for l=1:size(p,1)
                relK = -1:0.01:0;
                firstDeriv = halleyft(rho_d,polyval(p(l,:),relK));
                lastNaN = find(isnan(firstDeriv),1,'last');
                if(~isempty(lastNaN))
                    relK = relK(lastNaN:end);
                    secondDeriv = halleyft(rho_dd(:,lastNaN:end),firstDeriv([lastNaN+1 lastNaN+1:end]));
                    secondDerivValueAtFirst = interpft1([0 2*pi],rho_d(:,lastNaN:end),secondDeriv);
                    zeroCrossPos = find(diff(secondDerivValueAtFirst > 0),1,'first');
                    if(~isempty(zeroCrossPos))
                        zeroCross{r,c,K}(l,1) = K+relK(zeroCrossPos) - secondDerivValueAtFirst(zeroCrossPos)/diff(secondDerivValueAtFirst(zeroCrossPos+[0 1]))*0.01;
                        zeroCross{r,c,K}(l,2) = secondDeriv(zeroCrossPos) - secondDerivValueAtFirst(zeroCrossPos)/diff(secondDerivValueAtFirst(zeroCrossPos+[0 1]))*diff(secondDeriv(zeroCrossPos+[0 1]));
                        zeroCross{r,c,K}(l,2) = zeroCross{r,c,K}(l,2)/2;
                    end
                end
            end
        end
    end
end
save('zeroCross.mat','zeroCross');
        

%% Figure out spline

d = orientationSpace.diffusion.orientationMaximaDerivatives(rho,7:8,3,m(:,7:8)*2,pi*2);
md = cat(3,m(:,7:8)*2,d);
coefs = slvblk(colmat,joinColumns(squeeze(md(1,:,:)).')).';
sp = spmak(knots,coefs);


