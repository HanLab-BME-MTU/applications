function [handles,fhandle]=printKinPoleRef(kinTrack,EB3Tracks,varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addParameter('EB3ColorAfter',ones(1,length(EB3Tracks)));
ip.addParameter('valueRange',[]);
ip.addParameter('colormap',[], @islogical);
ip.parse(varargin{:});
p=ip.Results;



[handles,~,fhandle]=setupFigure(1,1,'AxesWidth',8,'AxesHeight',4,'DisplayMode', 'print');
hold(handles(1),'on');

scatter(handles,kinTrack.z,kinTrack.x,'r');
scatter(handles,0,0,'g');

if(~isempty(EB3Tracks))
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



for mIdx=1:length(EB3Tracks)
    mt=EB3Tracks(mIdx);

    %Project on the plan defined by the poleKin axis and the interpolar
    %axis.
    plot(handles,mt.z,mt.x,'Color',localCM(EB3ColorIdx(mIdx),:));
    ylim(handles,[-2000 2000]);
    xlabel(handles,'Pole-Kinetochore axis (nm)');
    ylabel(handles,'Normal plane (nm)');
end
    colorbar(handles,'ticks',linspace(0,1,5),'tickLabels',valueRange(round(linspace(1,end,5))));
end
% print([outputDirProj 'kin' num2str(kIdx,'%03d') '.png'],'-dpng');

% hold(handles(1),'off');
% hold(handles(2),'off');
% close(fhandle);
end
