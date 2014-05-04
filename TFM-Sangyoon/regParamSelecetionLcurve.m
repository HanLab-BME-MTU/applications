function [reg_corner,ireg_corner,kappa]=regParamSelecetionLcurve(rho,eta,lambda,manualSelection)%,dataPath)
% [reg_corner,ireg_corner,kappa]=l_curve_corner(rho,eta,reg_param)
% returns l curve corner estimated using a maximum curvature (kappa) estimation 
% in log-log space
% rho is the misfit and eta is the model norm or seminorm
%
% INPUT
%   rho       - misfit
%   eta       - model norm or seminorm
%   reg_param - the regularization parameter
%   manualSelection - true if you want to iterate selection process to get
%   better reg parameter. (default: false)
%
% OUTPUT
%   reg_corner  - the value of reg_param with maximum curvature
%   ireg_corner - the index of the value in reg_param with maximum curvature
%   kappa       - the curvature for each reg_param

if nargin<4
    manualSelection = false;
end

% transform rho and eta into log-log space
x=log(rho);
y=log(eta);
numCutPoints = 0;

% show the l curve and make sure this is fittable with 5th order polynomial
h=figure; set(h,'Position',[1000,100,400,1000])
subplot(3,1,1),plot(x,y,'k')
xlabel('Residual Norm ||Gm-d||_{2}');
ylabel('Simi-Norm ||Lm||_{2}');
x_cut = x((numCutPoints+1:end));
y_cut = y((numCutPoints+1:end));

kappa = diff(diff(y_cut)./diff(x_cut))./diff(x_cut(1:end-1));
subplot(3,1,2), plot(x_cut(1:end-2),kappa), title('curvature')
kappadiff = diff(kappa);
subplot(3,1,3), plot(x_cut(1:end-3),diff(kappa)),title('jerk')

% find a local maximum with three sections
nSections = 3;
nPoints = length(kappa);
p=0;
maxKappaCandIdx = [];
for ii=1:nSections
    [~, maxKappaIdx] = max(kappa(floor((ii-1)*nPoints/nSections)+1:floor(ii*nPoints/nSections))); % this is right at the L-corner which is usually over-smoothing
    maxKappaIdx = maxKappaIdx + floor((ii-1)*nPoints/nSections);
    % check if this is truly local maximum
    if maxKappaIdx>1 && maxKappaIdx<nPoints && (kappa(maxKappaIdx)>kappa(maxKappaIdx-1) && kappa(maxKappaIdx)>kappa(maxKappaIdx+1))
        p=p+1;
        maxKappaCandIdx(p) = maxKappaIdx;
    end
end
if length(maxKappaCandIdx)==1
    maxKappaIdx = maxKappaCandIdx(1);
elseif length(maxKappaCandIdx)>1
    [~,tempIndex] = max(kappa(maxKappaCandIdx));% use the first one %max(maxKappaCandIdx);
    maxKappaIdx = maxKappaCandIdx(tempIndex);
elseif isempty(maxKappaCandIdx)
    error('there is no local maximum in curvature in the input lambda range');
end
% [~, maxKappaDiffIdx] = max(kappadiff(1:maxKappaIdx)); %  this is steepest point right before L-corner. This is usually too small.
% find an index at kappa = 0 before maxKappaIdx
ireg_corner= numCutPoints+maxKappaIdx;%round((maxKappaIdx+maxKappaDiffIdx)/2); % thus we choose the mean of those two points.
reg_corner = lambda(ireg_corner);
disp(['Initial L-corner regularization parameter value: ' num2str(reg_corner) '.'])

% poly-fit version
p=polyfit(maxKappaIdx-5:maxKappaIdx+5,kappa(maxKappaIdx-5:maxKappaIdx+5)',2);
if p(1)<0
    ireg_corner = -p(2)/(2*p(1));
    q=polyfit(maxKappaIdx-5:maxKappaIdx+5,lambda(maxKappaIdx-5:maxKappaIdx+5),1);
    reg_corner = polyval(q,ireg_corner);
    disp(['Sub-knot resolution L-corner regularization parameter value: ' num2str(reg_corner) '.'])
else
    disp('The corner''''s L-corner does not have positive curvature')
    reg_corner = lambda(ireg_corner);
end

if manualSelection
    subplot(3,1,1), hold on,plot(x(ireg_corner),y(ireg_corner),'ro')
    text(x(ireg_corner),1.05*y(ireg_corner),...
        ['    ',num2str(lambda(ireg_corner),'%5.3e')]);

    subplot(3,1,2), hold on,plot(x_cut(ireg_corner),kappa(ireg_corner),'ro')
    subplot(3,1,3), hold on,plot(x_cut(ireg_corner),kappadiff(ireg_corner),'ro')

    poly5ivity = input('Is the curve going down with two concaveness (y/n)?','s');
    while poly5ivity == 'n'
        numCutPoints = input('how many entry points do you want to eliminate from the beginning?');
        subplot(3,1,1),plot(x(numCutPoints+1:end),y(numCutPoints+1:end),'k')
        x_cut = x((numCutPoints+1:end));
        y_cut = y((numCutPoints+1:end));

        kappa = diff(diff(y_cut)./diff(x_cut))./diff(x_cut(1:end-1));
        subplot(3,1,2), plot(x_cut(1:end-2),kappa)
        kappadiff = diff(kappa);
        subplot(3,1,3), plot(x_cut(1:end-3),diff(kappa))

        p=0;
        maxKappaCandIdx = [];
        nPoints = length(kappa);
        for ii=1:nSections
            [~, maxKappaIdx] = max(kappa(floor((ii-1)*nPoints/nSections)+1:floor(ii*nPoints/nSections))); % this is right at the L-corner which is usually over-smoothing
            maxKappaIdx = maxKappaIdx + floor((ii-1)*nPoints/nSections);
            % check if this is truly local maximum
            if maxKappaIdx>1 && maxKappaIdx<nPoints && (kappa(maxKappaIdx)>kappa(maxKappaIdx-1) && kappa(maxKappaIdx)>kappa(maxKappaIdx+1))
                p=p+1;
                maxKappaCandIdx(p) = maxKappaIdx;
            end
        end
        if length(maxKappaCandIdx)==1
            maxKappaIdx = maxKappaCandIdx(1);
        elseif length(maxKappaCandIdx)>1
            maxKappaIdx = max(maxKappaCandIdx);
        elseif isempty(maxKappaCandIdx)
            error('there is no local maximum in curvature in the input lambda range');
        end
    %     [~, maxKappaDiffIdx] = max(kappadiff(1:maxKappaIdx)); %  this is steepest point right before L-corner. This is usually too small.
        % find an index at kappa = 0 before maxKappaIdx
        ireg_corner= numCutPoints+maxKappaIdx;%round((maxKappaIdx+maxKappaDiffIdx)/2); % thus we choose the mean of those two points.
        subplot(3,1,1), hold on,plot(x(ireg_corner),y(ireg_corner),'ro')
        text(x(ireg_corner),1.1*y(ireg_corner), ['    ',num2str(lambda(ireg_corner),'%5.3e')]);

        subplot(3,1,2), hold on,plot(x_cut(ireg_corner-numCutPoints),kappa(ireg_corner-numCutPoints),'ro')
        subplot(3,1,3), hold on,plot(x_cut(ireg_corner-numCutPoints),kappadiff(ireg_corner-numCutPoints),'ro')

        poly5ivity = input('Is the curve going down with two concaveness (y/n)?','s');
    end
    reg_corner = lambda(ireg_corner);
end
% % fit it in 5th order polynomial - fitting with polynomial is dangerous!
% f = fit(x(numCutPoints+1:end), y(numCutPoints+1:end),  'poly5');
% hold on, plot(f)
% [~, fxx] = differentiate(f,x(numCutPoints+1:end));
% 
% % third derivative for finding local peak in curvature
% p = [f.p1 f.p2 f.p3 f.p4 f.p5 f.p6];
% 
% pxxx = polyder(polyder(polyder(p)));
% peaks2 = roots(pxxx);
% peak = peaks2(1);
% [~,peakIdx]=min(abs(x(numCutPoints+1:end)-peak)); % this is right at the L-corner which is over smoothing
% 
% pxxxx = polyder(polyder(polyder(polyder(p))));
% curvCurv = roots(pxxxx); % curvature of curvature
% [~,curvCurvIdx]=min(abs(x(numCutPoints+1:end)-curvCurv)); % this is steepest point right before L-corner. This is usually too small.
% 
% ireg_corner= round((peakIdx+curvCurvIdx)/2); % thus we choose the mean of those two points.
% reg_corner = lambda(numCutPoints+ireg_corner);
% kappa = fxx;
% 
% hold on
% plot(x(ireg_corner+numCutPoints),y(ireg_corner+numCutPoints),'ro')
% % print(h,[dataPath filesep 'Lcurve.eps'],'-depsc')
% close(h)
% find a positive peak in curvature (diff based, discrete)



