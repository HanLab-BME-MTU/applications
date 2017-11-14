function [hfig, B, R] = kymograph(R,K,r,c,out)
% Kymograph - Make kymograph figure for paper

INTERPOLATION_INTERVAL = 5;
PLOT_INTERVAL = 1;

if(nargin < 1)
    I = orientationSpace.paper.loadLaminDemoImage;
    F = OrientationSpaceFilter.constructByRadialOrder(1/2/pi./2,1,8,'none');
    R = F*I;
    r = 622;
    c = 364;
end

if(nargin < 2)
    K = 8:-0.1:1;
end


if(isnumeric(R))
    rho = R;
else
    if(nargin < 4 && ~exist('c','var'))
        r = round(size(R.a,1)/2);
        c = round(size(R.a,2)/2);
    end
    rho = R.getResponseAtOrderFTatPoint(r,c,K);
end

if(nargin < 5)
    %% Conditioning
    rhoh = fft(rho);
    rhoh(1,:) = 0;
    rhoh(abs(rhoh) < eps*1e3) = 0;
    rho = ifft(rhoh);
    out = interpft_extrema(rhoh,1,[],[],false);

    [out,out_events] = orientationSpace.diffusion.alignExtrema(out);
    out = orientationSpace.diffusion.unwrapExtrema(out,out_events);
end

%% Find end
interpolationIdx = 1:INTERPOLATION_INTERVAL:length(K);
lastK = cumsum(~isnan(out(:,interpolationIdx)),2);
% lastK = lastK(:,1:INTERPOLATION_INTERVAL:length(K));
lastK = max(lastK,[],2);
lastK = interpolationIdx(lastK).';
lastInd = sub2ind(size(out),(1:size(out,1)).',lastK);

newtonfig = figure;
[xg,Kg] = orientationSpace.diffusion.newtonBPprotoSimple(R,r,c);
close(newtonfig);

xgAligned = orientationSpace.diffusion.alignExtrema([out(lastInd) xg.']);

Kg_aligned = Kg;
for i = 1:length(Kg)
    Kg_aligned(i) = Kg(xgAligned(i,2) == xg.');
end

xgAligned(Kg_aligned < K(end),2) = xgAligned(Kg_aligned < K(end),1);
Kg_aligned(Kg_aligned < K(end)) = K(end);

%% Calculate derivatives

derivOrder = 3;
lm = out;
period = 2*pi;
[ dnm_dKn, dnm_dtn ] = orientationSpace.diffusion.orientationMaximaDerivatives( rho, K, derivOrder, lm , period);
dm_dK = dnm_dKn(:,:,1);
d2m_dK2 = dnm_dKn(:,:,2);

%% Calculate bifurcationDerivatives
interpolateBifurcation = Kg_aligned ~= K(end);
% xgAligned = xgAligned(interpolateBifurcation,:);
% Kg_aligned = Kg_aligned(interpolateBifurcation);
% lastK = lastK(interpolateBifurcation);

lm = xgAligned(:,1);
% dnK_dmn(1,:,:) are the 1st and 2nd derivatives at the last sampled K
[ dtn_dnm, maximaDerivatives, dnK_dmn ] = orientationSpace.diffusion.orientationMaximaTimeDerivatives( rho(:,lastK), K(lastK), derivOrder, lm.' , period);

% We might not need to calculate the derivatives at the bifurcation because
% we know them...
% dnK_dmn(2,:,1) == 0
% dnK_dmn(1,:,1) == -1/D * second derivative t to K factor
rho2 = R.getResponseAtOrderFTatPoint(r,c,Kg_aligned);
% Conditioning
rhoh2 = fft(rho2);
rhoh2(1,:) = 0;
rhoh2(abs(rhoh2) < eps*1e3) = 0;
rho2 = ifft(rhoh2);

lm = xgAligned(:,2);
% dnK_dmn(2,:,:) are the 1st and 2nd derivatives at the bifurcation K
[ dtn_dnm, maximaDerivatives, dnK_dmn(2,:,:) ] = orientationSpace.diffusion.orientationMaximaTimeDerivatives( rho2, Kg_aligned, derivOrder, lm.' , period);


%% Do interpolation

sp2c = cell(1,size(out,1));
xqc = cell(size(sp2c));

for trackNum = 1:size(out,1)
    try
%         figure;
%         title(['Track ' num2str(trackNum)]);
        idx_select = fliplr(1:INTERPOLATION_INTERVAL:length(K));
        track = out(trackNum,idx_select);
        idx_select = idx_select(~isnan(track));
        track = track(~isnan(track));

%         x = joinColumns(repmat(K(idx_select),2,1));


%         ygrad = joinColumns([track; -outg(trackNum,idx_select)]);
%         spgrad = spapi(augknt(K(idx_select),4,2),x.',ygrad.');

%         y = joinColumns([track; dm_dK(trackNum,idx_select)]);
%         sp = spapi(augknt(K(idx_select),4,2),x.',y.');

        x2 = joinColumns(repmat(K(idx_select),3,1));
        y2 = joinColumns([track; dm_dK(trackNum,idx_select); d2m_dK2(trackNum,idx_select)]);
        sp2 = spapi(optknt(x2.',5),x2.',y2.');
        sp2c{trackNum} = sp2;
%         sp2c{trackNum} = sp;


        xq = K(idx_select(1)):0.01:K(idx_select(end));
        xqc{trackNum} = xq;
%                 hold on;
%         plot(xq,spval(sp,xq))
%         plot(xq,spval(spgrad,xq))
%         plot(xq,spval(sp2,xq),'--')
%         plot(K,out(trackNum,:),'o');
%         plot(K(idx_select),out(trackNum,idx_select),'o','MarkerFaceColor','k');
%         keyboard
    catch err
    end
end

%% Do birfurcation interpolation

sppc = cell(1,size(out,1));
xThetaC = cell(1,size(out,1));
for trackNum = find(interpolateBifurcation)
    try
        if(xgAligned(trackNum,1) > xgAligned(trackNum,2))
            xTheta = joinColumns(repmat(fliplr(xgAligned(trackNum,:)),derivOrder+1,1));
            xThetaC{trackNum} = [xgAligned(trackNum,2):0.01:xgAligned(trackNum,1) xgAligned(trackNum,1)];
            yK = joinColumns(fliplr([K(lastK(trackNum)) Kg_aligned(trackNum); [squeeze(dnK_dmn(1,trackNum,:)) squeeze(dnK_dmn(2,trackNum,:))] ]));
            spp = spapi(optknt(xTheta.',5),xTheta.',yK.');
            sppc{trackNum} = spp;
        else
            xTheta = joinColumns(repmat(xgAligned(trackNum,:),derivOrder+1,1));
            xThetaC{trackNum} = [xgAligned(trackNum,1):0.01:xgAligned(trackNum,2) xgAligned(trackNum,2)];
            yK = joinColumns([K(lastK(trackNum)) Kg_aligned(trackNum); [squeeze(dnK_dmn(1,trackNum,:)) squeeze(dnK_dmn(2,trackNum,:))] ]);
            spp = spapi(optknt(xTheta.',5),xTheta.',yK.');
            sppc{trackNum} = spp;
        end
    catch err
    end
end

% orientationSpace.diffusion.diffusion_test

A = interpft(rho,360);
B = A;
for i=1:size(A,2)
     B(:,i) = mat2gray(A(:,i));
end

%% Plot

hfig = figure;
imagesc((0:359)/2,K,B.');
axis xy
colormap(bone);
hold on;
hp1 = plot(out(:,1:PLOT_INTERVAL:end)/pi/2*180,K(1:PLOT_INTERVAL:end),'.');
hp2 = plot(out(:,1:INTERPOLATION_INTERVAL:end)/pi/2*180,K(1:INTERPOLATION_INTERVAL:end),'o');
hp3 = gobjects(size(hp1));
for i=1:size(xgAligned,1)
    hp3(i) = plot(xgAligned(i,2)/pi/2*180,Kg_aligned(i),'s');
end
for i=1:length(hp1)
    try
        set(hp2(i),'Color',hp1(i).Color);
        set(hp3(i),'Color',hp1(i).Color);
    catch err
    end
end
for i=1:length(hp2)
    if(~isempty(sp2c{i}))
        plot(spval(sp2c{i},xqc{i})/2/pi*180,xqc{i},'Color',hp2(i).Color);
    end
    if(~isempty(sppc{i}))
        plot(xThetaC{i}/2/pi*180,spval(sppc{i},xThetaC{i}),'--','Color',hp2(i).Color);
    end
%     plot(spval(sp2c{i},xqc{i})/2/pi*180,xqc{i},'Color','m')
end
hcb = colorbar; 
hcb.TickLabels{end} = '1 Max';
hcb.TickLabels{1} = '0 Min';
hcb.Label.String = 'Relative Response to Max and Min per K';
ylabel('K');
xlabel('Orientation (degrees)');

end