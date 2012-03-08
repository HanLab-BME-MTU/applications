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
validFlowFrames = find(~cellfun(@isempty,flow));
assert(isequal(unique(diff(validFlowFrames)),1),['Flow must be measured on consecutive'...
    ' frames to be analyzed']);
startFrames = validFlowFrames(1:end-nAvg+1);
midFrames = startFrames+fix(nAvg/2);
preFrames=1:midFrames(1)-1;
postFrames = midFrames(end)+1:numel(flow)+1;

% Initialize empty output
Md=cell(1,numel(flow)+1);
Ms=cell(1,numel(flow)+1);
E=cell(1,numel(flow)+1);
S=cell(1,numel(flow)+1);

% Interpolate vector field
for i=startFrames
    Md{i+fix(nAvg/2)}=vectorFieldAdaptInterp(vertcat(flow{i:i+nAvg-1}),...
        flow{i+fix(nAvg/2)}(:,1:2),corrLength,[],'strain');
end
Md(preFrames)=Md(midFrames(1));
Md(postFrames)=Md(midFrames(end));

if ~ip.Results.noise, return; end

% Return cell array of noise vectors
for i=midFrames
    Ms{i}=horzcat(Md{i}(:,3:4),flow{i}(:,3:4));
end
Ms(preFrames)=Ms(midFrames(1));
Ms(postFrames)=Ms(midFrames(end));

if ~ip.Results.error, return; end

for i=midFrames
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

E(preFrames)=E(midFrames(1));
E(postFrames)=E(midFrames(end));

S(preFrames)=S(midFrames(1));
S(postFrames)=S(midFrames(end));