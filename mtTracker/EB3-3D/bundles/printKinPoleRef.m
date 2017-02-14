function [handles,fhandle]=printKinPoleRef(kinTrack,EB3Tracks)
% $$$ ip = inputParser;
% $$$ ip.CaseSensitive = false;
% $$$ ip.KeepUnmatched = true;
% $$$ ip.addRequired('MD',@(MD) isa(MD,'MovieData'));
% $$$ ip.addParameter('kinTracks',[],@(x) isa(x,'Tracks'));
% $$$ ip.addParameter('printAll',false, @islogical);
% $$$ ip.parse(kinTrack,EB3Tracks,varargin{:});
% $$$ p=ip.Results;


[handles,~,fhandle]=setupFigure(1,1,'AxesWidth',8,'AxesHeight',4,'DisplayMode', 'print');
hold(handles(1),'on');

for poleId=1:2
    scatter(handles,kinTrack.z,kinTrack.x,'r');
    scatter(handles,0,0,'g');
end

for mIdx=1:length(EB3Tracks)
    mt=EB3Tracks(mIdx);
    
    %Project on the plan defined by the poleKin axis and the interpolar
    %axis.
    plot(handles,mt.z,mt.x,'b-');
    ylim(handles,[-2000 2000])
    xlabel(handles,'Pole-Kinetochore axis (nm)')
    ylabel(handles,'Normal plane (nm)')
end

% print([outputDirProj 'kin' num2str(kIdx,'%03d') '.png'],'-dpng');

% hold(handles(1),'off');
% hold(handles(2),'off');
% close(fhandle);
end
