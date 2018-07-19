function [handles]=printMTTipsPoleRef(kinTrack,EB3Tracks,MTTipsIdx,varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addParameter('EB3ColorAfter',ones(1,length(EB3Tracks)));
ip.addParameter('valueRange',[]);
ip.addParameter('normAxis',false);
ip.addParameter('handle',[]);

ip.addParameter('colormap',[], @islogical);
ip.parse(varargin{:});
p=ip.Results;

if(isempty(EB3Tracks))
    p.EB3ColorAfter=1;
else  

minColoredValue=min(p.EB3ColorAfter);
maxColoredValue=max(p.EB3ColorAfter);

valueRange=p.valueRange;
if(isempty(p.valueRange))
    valueRange=floor(minColoredValue):ceil(maxColoredValue);
else
    valueRange=valueRange(1):valueRange(2);
end

[EB3ColorIdx] = discretize(p.EB3ColorAfter,[min(minColoredValue,min(valueRange)) valueRange  max(maxColoredValue,max(valueRange))+1]);

localCM=p.colormap;
if(isempty(p.colormap))
  localCM=jet(max(EB3ColorIdx));
end

handles=p.handle;
if(isempty(p.handle))
[handles,~,fhandle]=setupFigure(1,1,'AxesWidth',8,'AxesHeight',4,'DisplayMode', 'print');
end
hold(handles(1),'on');
if(~p.normAxis)
scatter(handles,kinTrack.z,kinTrack.x,'r');
else
    scatter(handles,ones(size(kinTrack.z)),kinTrack.x,'r');
end
scatter(handles,0,0,'g');

X=nan(1,length(EB3Tracks));
Y=nan(1,length(EB3Tracks));
for mIdx=1:length(EB3Tracks)
    mt=EB3Tracks(mIdx);
    if(~p.normAxis)
        X(mIdx)=mt.z(MTTipsIdx(mIdx));
    else
        X(mIdx)=mt.z(MTTipsIdx(mIdx))/kinTrack.z(mt.f(MTTipsIdx(mIdx))==kinTrack.f);
    end
    Y(mIdx)=mt.x(MTTipsIdx(mIdx));
    %Project on the plan defined by the poleKin axis and the interpolar
    %axis.
end


scatter(handles,X,Y,10,localCM(EB3ColorIdx,:));
if(p.normAxis)
    xlim(handles,[0 1.2]);
end
ylim(handles,[-2000 2000]);
xlabel(handles,'Pole-Kinetochore axis (nm)');
ylabel(handles,'Normal plane (nm)');
colorbar(handles,'ticks',linspace(0,1,5),'tickLabels',valueRange(round(linspace(1,end,5))));

% print([outputDirProj 'kin' num2str(kIdx,'%03d') '.png'],'-dpng');

% hold(handles(1),'off');
% hold(handles(2),'off');
% close(fhandle);
end
