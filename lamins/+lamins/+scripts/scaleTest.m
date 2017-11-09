% FF = OrientationSpaceRidgeFilter(1./2/pi./3.1,[],8,'energy');
% FFR = FF*delta;
% template = FFR.a(:,:,1);
% template(template < 0) = 0;
% template = repmat(normpdf(-49:50,0,2)+5,100,1);
template = imgaussfilt(delta,pi,'FilterSize',61);
% template = bsxfun(@times,normpdf(-49:50,0,1.89),normpdf(-49:50,0,5)');
start = 1;
finish = 8;
npts = 513;
TOL = 1e-6;
s = chebpts(npts,[start finish]);
iters = log2(length(s)-1)-1;
F2 = OrientationSpaceRidgeFilter(1./2/pi./s,[],0,'none');
F2([1 npts]).setupFilter(size(template));
F2R([1 npts],1) = F2([1 npts])*template;
% F2([1 npts]).clearTransients;
clear test
figure;
for i=1:iters+1
% for p=2.^(iters:-1:0)
    p = 2.^(iters+1-i);
    idx = (1+p):p*2:npts;
    F2(idx).setupFilter(size(template));
    F2R(idx) = F2(idx)*template;
%     F2(idx).clearTransients;
    idx = 1:p:npts;
%     F2R = F2*template;
    % N = sqrt(sum(abs(F2(1).F(:)).^2));
    % for i=1:17; F2(i).F = F2(i).F / N; end;
    % F2R = F2*FFR.a(:,:,1);
    % F2Ras = F2R.getArraySpace;
    % cf = chebfun(squeeze(F2Ras(50,50,1,:)),[1 4]);
    % plot(cf)
    % hold on;

    % for i=1:length(F2); F2(i).F = F2(i).F / sqrt(sum(abs(F2(i).F(:)).^2)); end;
%     for i=1:length(F2); F2(i).F = F2(i).F .* s(i); end;
%     F2R = F2*FFR.a(:,:,1);
    F2Ras = F2R(idx).getArraySpace;
    r = squeeze(F2Ras(50,50,1,:));
    r = r./real(F2(idx).getEnergy).^0.6105;
    cf = chebfun(r,[start finish]);
    plot(cf,'-o','MarkerSize',i*3+3)
    hold on
%     plot(cf,'o','MarkerSize',i*3+5)
    [~,test(i)] = max(cf)
    plot([test(i) test(i)],ylim);
    pause(0.1);
    if(i > 1 && abs(test(i) - test(i-1)) < TOL)
        break;
    end
    pause(1);
end

clear F2R

figure; plot(test);


% F = OrientationSpaceRidgeFilter(1./2/pi./3.1428,[],8); F.setupFilter([ 100 100]); F.F = F.F*3.1428;
% FR = F*template;
% FR.a(50,50,1)
