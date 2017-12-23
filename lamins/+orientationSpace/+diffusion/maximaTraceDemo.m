function [q] = maximaTraceDemo()

gcp;

delta = 0.5;
K = 8:-delta:3;
threshold = pi/60;

I = orientationSpace.paper.loadLaminDemoImage;
F = OrientationSpaceFilter.constructByRadialOrder(1/2/pi./2,1,8,'none');
R = F*I;
maxima = R.getRidgeOrientationLocalMaxima;
disp('Doing Maxima Trace');
maximaTrace = R.traceAllLocalMaxima(K,maxima);

disp('Calculating Derivatives');
% res = cell(1,length(K));
% for ki=1:length(K)
%   res{ki} = real(shiftdim(R.getResponseAtOrderFT(K(ki),2).a,2));
% end
tic
parfor (ki=1:length(K),6)
    res = real(shiftdim(R.getResponseAtOrderFT(K(ki),2).a,2));
    lm = shiftdim(maximaTrace(:,:,:,ki),2)*2;
    lmd = orientationSpace.diffusion.orientationMaximaDerivatives(res,K(ki),2,lm);
    maximaTraceDeriv(:,:,:,ki) = shiftdim(lmd(:,:,:,1),1);
    maximaTraceDeriv2(:,:,:,ki) = shiftdim(lmd(:,:,:,2),1);
end
toc

disp('Interpolating derivatives and finding threshold');
tic
maximaTraceDerivIsNaN = isnan(maximaTraceDeriv);
maximaTraceDerivNotNaN = ~maximaTraceDerivIsNaN;
maximaTraceDerivIndexed = cumsum(maximaTraceDerivNotNaN,4);
maximaTraceDerivCount = sum(maximaTraceDerivNotNaN,4);
maximaTraceDerivCountBelowThreshold = sum(maximaTraceDerivNotNaN & abs(maximaTraceDeriv) <= threshold,4);
maximaTraceDerivCountAboveThreshold = maximaTraceDerivCount - maximaTraceDerivCountBelowThreshold;
doInterpolation = maximaTraceDerivCountBelowThreshold >  0 & maximaTraceDerivCountAboveThreshold  > 0;
findBifurcation = maximaTraceDerivCountBelowThreshold >  0 & maximaTraceDerivCountAboveThreshold == 0;
acceptCurrent   = maximaTraceDerivCountBelowThreshold == 0 & maximaTraceDerivCountAboveThreshold  > 0;
maximaTraceDerivIndexed(maximaTraceDerivIsNaN) = NaN;
maximaTraceDerivIndexed(abs(maximaTraceDeriv) > threshold) = NaN;
ind = max(maximaTraceDerivIndexed,[],4);
indIsNaN = isnan(ind);
indIsNotNaN = ~indIsNaN;
indNoNaN = ind;
indNoNaN(indIsNaN) = 1;

% [~,ind] = max(cumsum(abs(maximaTraceDeriv) < threshold,4),[],4);
[q.X,q.Y,q.Z] = meshgrid(1:1024,1:1024,1:7);
q.maximaTraceSelect = maximaTrace(sub2ind(size(maximaTrace),q.Y,q.X,q.Z,indNoNaN));
R3 = R.getResponseAtOrderFT(3,2);
q.nlms = R3.nonLocalMaximaSuppressionPrecise(q.maximaTraceSelect);
q.maximaTrace = maximaTrace;
q.maximaTraceDeriv = maximaTraceDeriv;
q.maximaTraceDeriv2 = maximaTraceDeriv2;

q.maximaTraceData = cat(5,q.maximaTrace*2,q.maximaTraceDeriv,q.maximaTraceDeriv2);
q.ind = ind;

k = 4;
tau = [zeros(1,k-1) delta*ones(1,k-1)];
knots = optknt(tau,k);
blockmat = spcol(knots,k,tau,'sl');
r = 622;
c = 364;

% ind = ind(indIsNotNaN);
% sMaxed = ind == maximaTraceDerivCount;
% s = ~sMaxed & ~isnan(ind);
% s = true(size(ind));
s = indIsNotNaN;
ind = ind(s);
selectPrev = ind == maximaTraceDerivCount(s) & ind > 1;
selectNext = ~selectPrev;
ind(:,2) = NaN;
ind(selectPrev,:) = bsxfun(@plus,ind(selectPrev,1),[ 0 -1]);
ind(selectNext,:) = bsxfun(@plus,ind(selectNext,1),[ 1 0]);
sn = sum(s(:));

sz = size(q.maximaTraceData);
q.maximaTraceData = reshape(q.maximaTraceData,prod(sz(1:3)),sz(4),sz(5));
q.maximaTraceData = q.maximaTraceData(s,:,:);

linear_ind = sub2ind(size(q.maximaTraceData),repmat((1:sn).',1,2,3),repmat(ind,1,1,3),repmat(shiftdim(1:3,-1),sn,2,1));
linear_ind = permute(linear_ind,[3 2 1]);
linear_ind = reshape(linear_ind,6,size(linear_ind,3));

intervalData = q.maximaTraceData(linear_ind);
coefs = slvblk(blockmat,intervalData).';

sp = spmak(knots,coefs);
spd = fnder(sp);
% spd.dim = 1;
idxPos = spd.coefs(:,1) >  threshold & spd.coefs(:,5) <  threshold;
idxNeg = spd.coefs(:,1) < -threshold & spd.coefs(:,5) > -threshold;
spd.coefs(idxPos,:) = spd.coefs(idxPos,:) - threshold;
spd.coefs(idxNeg,:) = spd.coefs(idxNeg,:) + threshold;

pp = sp2pp(spd);
discriminant = pp.coefs(:,2).^2-4*pp.coefs(:,1).*pp.coefs(:,3);
dp = discriminant > 0;
sqrt_discriminant = sqrt(discriminant(dp));
pproots = NaN(size(pp.coefs,1),2);
pproots(dp,1) = (-pp.coefs(dp,2)+sqrt_discriminant)/2./pp.coefs(dp,1);
pproots(dp,2) = (-pp.coefs(dp,2)-sqrt_discriminant)/2./pp.coefs(dp,1);
pproots = reshape(pproots,size(pp.coefs,1)/3,3,2);
pproots(pproots > diff(pp.breaks)) = NaN;
pproots(pproots < 0) = NaN;
pproots = pproots + pp.breaks(1:3);
pproots = reshape(pproots,size(pproots,1),6);
pproots = max(pproots,[],2);

rr = NaN(sz(1:3));
rr(s) = pproots + K(ind(:,1)).';
q.rr = rr;
toc
%15.239661 seconds

% TODO: 

% coefsd = num2cell(spd.coefs,2).';
% % manipulate coefsd by threshold somehow
% spdc = struct2cell(spd);
% spdc = repmat(spdc,1,size(spd.coefs,1));
% spdc(3,:) = coefsd;
% spd = cell2struct(spdc,fields(spd));


% q.maximaTraceData = shiftdim(q.maximaTraceData,2);
% si = find(s);
% ii = sub2ind([size(q.maximaTraceData,2) size(q.maximaTraceData,3)],si,ind(:,2));
% yy = q.maximaTraceData(:,ii);
% ii = sub2ind([size(q.maximaTraceData,2) size(q.maximaTraceData,3)],si,ind(:,1));
% yy = [yy; q.maximaTraceData(:,ii)];

%% Code for analyzing spline at r,c
% 
% y = [q.maximaTrace(r,c,1,10)*2 q.maximaTraceDeriv(r,c,1,10) q.maximaTraceDeriv2(r,c,1,10) q.maximaTrace(r,c,1,9)*2 q.maximaTraceDeriv(r,c,1,9) q.maximaTraceDeriv2(r,c,1,9)].';
% coefs = slvblk(blockmat,y).';
% sp = spmak(knots,coefs);
% 
% 
% % sp = spapi(aptknt([0 0 0 0.5 0.5 0.5],4),[0 0 0 0.5 0.5 0.5],[q.maximaTrace(r,c,1,10)*2 q.maximaTraceDeriv(r,c,1,10) q.maximaTraceDeriv2(r,c,1,10) q.maximaTrace(r,c,1,9)*2 q.maximaTraceDeriv(r,c,1,9) q.maximaTraceDeriv2(r,c,1,9)]);
% pp = sp2pp(fnder(sp));
% pp.coefs(:,3) = pp.coefs(:,3)+pi/60;
% z = fnzeros(pp);

% coefs = slvblk(blockmat,yy).';
% z = zeros(2,size(coefs,1));
% for i=1:size(coefs,1)
%     sp = spmak(knots,coefs(i,:));
%     pp = sp2pp(fnder(sp));
%      pp.coefs(:,3) = pp.coefs(:,3)+threshold;
%      z(:,i) = fnzeros(pp);
% end

end