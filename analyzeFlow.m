function [Md,Ms,E,S] = analyzeFlow(flow,nAvg,corrLength,varargin)
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
% OUTPUT        Md          : interpolated vectors
%               Ms          : noise vectors
%               E           : cell array of errors between raw and
%               interpolated vector fields
%               S           : cell array of vector field signal to noise
%               ratios
%

% Aaron Ponti, May 15th, 2003
% Sebastien Besson, June 2011 (last modified Nov 2011)
% Adapted from fsmVectorAnalysis

% Check input
ip = inputParser;
ip.addRequired('flow',@iscell);
ip.addRequired('nAvg',@(x) isscalar(x) & mod(x,2)~=0);
ip.addRequired('corrLength',@isscalar);
ip.addParamValue('noise',1,@isscalar);
ip.addParamValue('error',1,@isscalar);
ip.parse(flow,nAvg,corrLength,varargin{:});

nImages= numel(flow)-nAvg+1;
frames = find(~arrayfun(@isempty,flow(1:nImages)));

% Initialize empty output
Md=cell(1,nImages);
Ms=cell(1,nImages);
E=cell(1,nImages);
S=cell(1,nImages);

% Interpolate vector field
Md(frames)=arrayfun(@(i) vectorFieldAdaptInterp(vertcat(flow{i:i+nAvg-1}),...
    flow{i+fix(nAvg/2)}(:,1:2),corrLength,[],'strain'),frames,...
    'UniformOutput',false);

if ~ip.Results.noise, return; end

% Return cell array of noise vectors
Ms(frames)=arrayfun(@(i) horzcat(Md{i}(:,3:4),flow{i+fix(nAvg/2)}(:,3:4)),frames,...
    'UniformOutput',false);

if ~ip.Results.error, return; end

for i=frames
    % Extract vectors
    d=Md{i}(:,3:4)-Md{i}(:,1:2); % Interpolated vector field
    v=flow{i+fix(nAvg/2)}(:,3:4)-flow{i+fix(nAvg/2)}(:,1:2); % Raw vector field
    s=Ms{i}(:,3:4)-Ms{i}(:,1:2); % Noise field
    
    % Vector lengths
    ld=sqrt(d(:,1).^2+d(:,2).^2);
    lv=sqrt(v(:,1).^2+v(:,2).^2);
    ls=sqrt(s(:,1).^2+s(:,2).^2);
    
    % Measure relative error on lengths between raw and interploated fields
    eL=abs(lv-ld)./ld;
    % Measure signal to noise ratio between interpolated and noise vectors
    SNR=ld./ls;
    
    % Angles - vectorize
    dy=d(:,1); dx=d(:,2); vy=v(:,1); vx=v(:,2);
    num=sqrt(vy.^2+vx.^2).*sqrt(dy.^2+dx.^2);
    num(num==0)=eps;
    eA=acos((vy.*dy+vx.*dx)./num)/pi;
    
    % Total error
    eT=sqrt(eL.^2+eA.^2);
    
    % Construct matrices for display
    E{i}=[flow{i+fix(nAvg/2)}(:,1:2) eT];
    S{i}=[flow{i+fix(nAvg/2)}(:,1:2) SNR];
end
        
%     
%     if ERROR_CALC==1
% 
%         fprintf(1,'File name                        : %s\n',imageFileList{c1});
%         fprintf(1,'Correlation length corrLength            : %d\n',corrLength_init);
%         fprintf(1,'Number of RAW vectors            : %d\n',size(flow,1));
%         fprintf(1,'Mean RAW vector length           : %2.4f +/- %2.4f (+/- %2.2f%%)\n',mean(lv),std(lv),100*std(lv)/mean(lv));
%         fprintf(1,'Median RAW vector length         : %2.4f\n',median(lv));
%         fprintf(1,'Mean INTERPOLATED vector length  : %2.4f +/- %2.4f (+/- %2.2f%%)\n',mean(ld),std(ld),100*std(ld)/mean(ld));
%         fprintf(1,'Median INTERPOLATED vector length: %2.4f\n',median(ld));
%         fprintf(1,'Mean NOISE vector length         : %2.4f +/- %2.4f (+/- %2.2f%%)\n',mean(ls),std(ls),100*std(ls)/mean(ls));
%         fprintf(1,'Median NOISE vector length       : %2.4f\n',median(ls));
%         snr=ld./ls;
%         fprintf(1,'Mean / median SNR                : %2.4f +/- %2.4f (+/- %2.2f%%) / %2.4f\n',mean(snr),std(snr),100*std(snr)/mean(snr),median(snr));
%         [n,h]=hist(snr,[0.5:max(snr)+0.5]);
%         n=cumsum(n); n=n/max(n);
%         indx=find(n>=0.95);indx=indx(1);
%         fprintf(1,'95%% of all vectors have SNR <= %d\n',indx);
%     end
    
