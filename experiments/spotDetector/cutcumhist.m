function cutpts=cutcumhist(mnp,ti,dataProperties)
%CUTCUMHIST separates the relevant spots in the cumulative histogram from false spots
%
% SYNOPSIS cutpts=cutcumhist(mnp,ti,plt)
%
% INPUT mnp : 'spotiness' list  = curvature*mean(intensity) (matrix #spts(t) x timepoints)
%            ti      : current time point (for mnp )  (scalar)
%            plt    : (Boolean) true: plot histogram , false: don't plot
% OUTPUT cutpts : index in mnp fto relevant spots

% c: 5/7/01	dT

%Define constants
%global CH_MAXNUMINTERV CH_MAXSLOPE;
CH_MAXNUMINTERV = dataProperties.CH_MAXNUMINTERV;
CH_MAXSLOPE = dataProperties.CH_MAXSLOPE;

%ipos = find(mnp(:,ti)>0);
%[mnp mnpidx]=sort(mnp(ipos));
[mnp mnpidx]=sort(nonzeros(mnp(:,ti)));
% calc sampling intervall
sint=sort(diff(mnp));
samplingIntr=sint(min(find(sint>0)));

%Check for path ological cases (only 1-4 mnp's)

% check for too small sampling intreval
if (mnp(end)-mnp(1))/samplingIntr>CH_MAXNUMINTERV
    samplingIntr=(mnp(end)-mnp(1))/CH_MAXNUMINTERV;
end;

% interpolate cumhist
cutpts=0;
x=mnp(1):samplingIntr:mnp(end);
uniqmnp=mnp(find(diff(mnp)));
uniqmnp(end+1)=mnp(end);
y=interp1(uniqmnp,1:length(uniqmnp),x);

%[MUHAT,SIGMAHAT,MUCI,SIGMACI] = NORMFIT(uniqmnp);
%y=normcdf(x,MUHAT,SIGMAHAT);

%find horizontal parts
dfy=diff(y)/samplingIntr;
[m i]=max(dfy); %get max histo-slope (i: index)
%find 0- MAXSPOTS
idx=[];
slope=CH_MAXSLOPE;
%slope of cumhist < MAX_SLOPE
zfy=find(dfy<slope); %cut histogram
idx =min(nonzeros(zfy.*(zfy>i))); 

%    idx=max(find(dfy>slope));

if(~isempty(idx))
    if idx <(length(x)-2) & length(mnp) > 10
        idx=idx+2;
    end;
    %cutpts=find(mnp>x(idx));
    cutpts=find(mnp>=max(x(idx),100));  % spotiness at least ==100
    cutpts=mnpidx(cutpts);
else
    cutpts=0;
end;

%plot preview panel if exist
if isfield(dataProperties,'previewHandle')
    axesH=dataProperties.previewHandle;
    axes(axesH);
    %reduce to the positive vals:
    ipos=find(x>=0);
    %plot(x(ipos),y(ipos),'b-');
    hold on;
    ipos=find(mnp>=0);
    plot(mnp(ipos),1:length(mnp(ipos)),'g+');
    ipos2 = find(mnp>=100);
    plot(mnp(ipos2),length(mnp(ipos))-length(mnp(ipos2))+1:length(mnp(ipos)),'b+');
    if ~isempty(idx)
        tp=find(mnp>=max(x(idx),100));
        plot(mnp(tp),length(mnp(ipos))-length(tp)+1:length(mnp(ipos)),'r+');
    end
end;