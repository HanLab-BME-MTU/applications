function [Md,Ms,E,S,stats] = analyzeFlow(flow,nAvg,corrLength,varargin)
% analyzeFlow reports several statistics on vector fields (obtained from FSM)
%
% SYNOPSIS      analyzeFlow(flow,nAvg,corrLength)
%               analyzeFlow(flow,nAvg,corrLength,'interpolate',1)
%               analyzeFlow(flow,nAvg,corrLength,'interpolate',1,'noise',1)
%               analyzeFlow(flow,nAvg,corrLength,'interpolate',1,'noise',1,'error',1)
%
%
% INPUT         flow        : cell array of Nx4 arrays containing the flow
%               nAvg        : number of frames for time averaging (must be odd)
%               corrLength  : correlation length of the interpolator
%                
% OUTPUT        Md      : interpolated vectors
%               Ms      : noise vectors
%               E       : cell array of errors between raw and interpolated 
%                         vector fields
%               S       : cell array of vector field signal to noise ratios
%               stats   : structure containing flow statistics
%

% Aaron Ponti, May 15th, 2003
% Sebastien Besson, June 2011 (last modified Nov 2011)
% Adapted from fsmVectorAnalysis

% Check input
ip = inputParser;
ip.addRequired('flow',@iscell);
ip.addRequired('nAvg',@(x) isscalar(x) & mod(x,2)~=0);
ip.addRequired('corrLength',@isscalar);
ip.addParamValue('noise',true,@isscalar);
ip.addParamValue('error',true,@isscalar);
ip.parse(flow,nAvg,corrLength,varargin{:});

% Get frame for interpolation
nInterpolations= numel(flow)-nAvg+1;
validFrames = find(~arrayfun(@isempty,flow(1:nInterpolations)));
startFrames=1:fix(nAvg/2);
endFrames = validFrames(end)+fix(nAvg/2)+1:numel(flow)+1;

% Initialize empty output
Md=cell(1,numel(flow)+1);
Ms=cell(1,numel(flow)+1);
E=cell(1,numel(flow)+1);
S=cell(1,numel(flow)+1);

% Interpolate vector field
for i=validFrames
    Md{i+fix(nAvg/2)}=vectorFieldAdaptInterp(vertcat(flow{i:i+nAvg-1}),...
        flow{i+fix(nAvg/2)}(:,1:2),corrLength,[],'strain');
end
Md(startFrames)=Md(fix(nAvg/2)+1);
Md(endFrames)=Md(nInterpolations+fix(nAvg/2));

if ~ip.Results.noise, return; end

% Return cell array of noise vectors
for i=validFrames+fix(nAvg/2)
    Ms{i}=horzcat(Md{i}(:,3:4),flow{i}(:,3:4));
end
Ms(startFrames)=Ms(fix(nAvg/2)+1);
Ms(endFrames)=Ms(nInterpolations+fix(nAvg/2));

if ~ip.Results.error, return; end

for i=validFrames+fix(nAvg/2)
    % Extract vectors
    v=flow{i}(:,3:4)-flow{i}(:,1:2); % Raw vector field
    d=Md{i}(:,3:4)-Md{i}(:,1:2); % Interpolated vector field    
    s=Ms{i}(:,3:4)-Ms{i}(:,1:2); % Noise field
    
    % Vector lengths
    ld=sqrt(d(:,1).^2+d(:,2).^2);
    lv=sqrt(v(:,1).^2+v(:,2).^2);
    ls=sqrt(s(:,1).^2+s(:,2).^2);
    
    % Measure relative error on lengths between raw and interploated fields
    eL=abs(lv-ld)./ld;
    % Measure signal to noise ratio between interpolated and noise vectors
    snr=ld./ls;
    
    % Angles - vectorize
    dy=d(:,1); dx=d(:,2); vy=v(:,1); vx=v(:,2);
    num=sqrt(vy.^2+vx.^2).*sqrt(dy.^2+dx.^2);
    num(num==0)=eps;
    eA=acos((vy.*dy+vx.*dx)./num)/pi;
    
    % Total error
    eT=sqrt(eL.^2+eA.^2);
    
    % Construct matrices for display
    E{i}=[flow{i}(:,1:2) eT];
    S{i}=[flow{i}(:,1:2) snr];
    
    stats.lv{i}=lv;
    stats.ld{i}=ld;
    stats.ls{i}=ls;
    stats.snr{i}=snr;
end

E(startFrames)=E(fix(nAvg/2)+1);
E(endFrames)=E(nInterpolations+fix(nAvg/2));
S(startFrames)=S(fix(nAvg/2)+1);
S(endFrames)=S(nInterpolations+fix(nAvg/2));


% fprintf(1,'Number of RAW vectors            : %d\n',size(flow,1));
% fprintf(1,'Mean RAW vector length           : %2.4f +/- %2.4f (+/- %2.2f%%)\n',mean(lv),std(lv),100*std(lv)/mean(lv));
% fprintf(1,'Median RAW vector length         : %2.4f\n',median(lv));
% fprintf(1,'Mean INTERPOLATED vector length  : %2.4f +/- %2.4f (+/- %2.2f%%)\n',mean(ld),std(ld),100*std(ld)/mean(ld));
% fprintf(1,'Median INTERPOLATED vector length: %2.4f\n',median(ld));
% fprintf(1,'Mean NOISE vector length         : %2.4f +/- %2.4f (+/- %2.2f%%)\n',mean(ls),std(ls),100*std(ls)/mean(ls));
% fprintf(1,'Median NOISE vector length       : %2.4f\n',median(ls));
% snr=ld./ls;
% fprintf(1,'Mean / median SNR                : %2.4f +/- %2.4f (+/- %2.2f%%) / %2.4f\n',mean(snr),std(snr),100*std(snr)/mean(snr),median(snr));
% [n,h]=hist(snr,0.5:max(snr)+0.5);
% n=cumsum(n); n=n/max(n);
% indx=find(n>=0.95);indx=indx(1);
% fprintf(1,'95%% of all vectors have SNR <= %d\n',indx);
% 
%     
