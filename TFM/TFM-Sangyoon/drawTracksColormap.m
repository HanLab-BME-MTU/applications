function h = drawTracksColormap(tracksNA,iFrame,Property,PropRange,Colormap)
% h = drawTracksColormap(tracksNA,iFrame,Property,PropRange,Colormap) draws
% trajectories of adhesions on top of an existing figure, colorcoded w.r.t.
% Property with PropRange.
% Usage:
%   
%   figure, imshow(


markerSize = 7;
circleMarkerSize=4;
lineWidth = 1;

% Index vectors per magnitude
nColors = size(Colormap,1);

intensity= arrayfun(@(x) x.(Property),tracksNA); %(u.^2+v.^2).^(1/2);
vColor=floor(scaleContrast(intensity,PropRange,[1 nColors]));
vColor(vColor<1)=1;
vColor(vColor>nColors)=nColors;
vIndex= unique(vColor);
hold on
% Create array of quiverplots
for i=1:numel(vIndex)
    idx = find(vColor==vIndex(i));

    xmat = cell2mat(arrayfun(@(x) x.xCoord(1:iFrame),tracksNA(idx),'UniformOutput',false));
    ymat = cell2mat(arrayfun(@(x) x.yCoord(1:iFrame),tracksNA(idx),'UniformOutput',false));
    if size(xmat,2)==1
        h{i}=plot(xmat',ymat','.','Color',Colormap(vIndex(i),:),'MarkerSize',markerSize);
    else
        h{i}=plot(xmat',ymat','Color',Colormap(vIndex(i),:),'LineWidth',lineWidth);
        plot(xmat(:,end)',ymat(:,end)','.','Color',Colormap(vIndex(i),:),'LineWidth',lineWidth/2,'MarkerSize',markerSize);
%         idAdhLogic = arrayfun(@(x) ~isempty(x.adhBoundary),tracksNA) & vColor==vIndex(i);
%         idAdhCur = arrayfun(@(x) ~isempty(x.adhBoundary{CurrentFrame}),tracksNA(idAdhLogic));
%         idAdh = find(idAdhLogic);
%         idAdhCur = idAdh(idAdhCur);
%         arrayfun(@(x) plot(x.adhBoundary{iFrame}(:,1),x.adhBoundary{CurrentFrame}(:,2), ...
%             'Color',Colormap(vIndex(i),:), 'LineWidth', 0.5),tracksNA(idAdhCur))
    end
end
