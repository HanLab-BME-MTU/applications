det=detLabRef(1)
[handles,~,fhandle]=setupFigure(2,2,3,'AxesWidth',6,'AxesHeight',6);
scatter(handles(1),det.xCoord(:,1),det.yCoord(:,1))
scatter(handles(2),det.xCoord(:,1),det.zCoord(:,1))
scatter(handles(3),det.yCoord(:,1),det.zCoord(:,1))

det=detLabRef(1).pole1
[handles,~,fhandle]=setupFigure(2,2,3,'AxesWidth',6,'AxesHeight',6);
scatter(handles(1),det.xCoord(:,1),det.yCoord(:,1))
scatter(handles(2),det.xCoord(:,1),det.zCoord(:,1))
scatter(handles(3),det.yCoord(:,1),det.zCoord(:,1))

det=detLabRef(1).pole2
[handles,~,fhandle]=setupFigure(2,2,3,'AxesWidth',6,'AxesHeight',6);
scatter(handles(1),det.xCoord(:,1),det.yCoord(:,1))
scatter(handles(2),det.xCoord(:,1),det.zCoord(:,1))
scatter(handles(3),det.yCoord(:,1),det.zCoord(:,1))

%%
detKins=kinTracksInliers(1).associatedDetectP1;
X=[];Y=[];Z=[];
for i=1:length(detKins)
    if(~isempty(detKins(i)))
        X=[X;detKins(i).xCoord];
        Y=[Y;detKins(i).yCoord];
        Z=[Z;detKins(i).zCoord];
    end
end
[handles,~,fhandle]=setupFigure(2,2,3,'AxesWidth',6,'AxesHeight',6);
scatter(handles(1),X(:,1),Y(:,1));
scatter(handles(2),X(:,1),Z(:,1));
scatter(handles(3),Y(:,1),Z(:,1));
xlim(handles(1),[0,MD.pixelSize_*MD.imSize_(1)]);
ylim(handles(1),[0,MD.pixelSize_*MD.imSize_(2)]);

xlim(handles(2),[0,MD.pixelSize_*MD.imSize_(1)]);
ylim(handles(2),[0,MD.pixelSizeZ_*MD.zSize_]);

xlim(handles(3),[0,MD.pixelSize_*MD.imSize_(2)]);
ylim(handles(3),[0,MD.pixelSizeZ_*MD.zSize_]);

%%

