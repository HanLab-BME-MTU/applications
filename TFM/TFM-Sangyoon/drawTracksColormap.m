function [h,PropRange] = drawTracksColormap(tracksNA,iFrame,Property,PropRange,Colormap,adhAnalProc)
% h = drawTracksColormap(tracksNA,iFrame,Property,PropRange,Colormap) draws
% trajectories of adhesions on top of an existing figure, colorcoded w.r.t.
% Property with PropRange.
% Usage:
%   
%   figure, imshow(


markerSize = 7;
circleMarkerSize=4;
lineWidth = 1;
p=adhAnalProc.funParams_;

% Index vectors per magnitude
nColors = size(Colormap,1);

intensity= arrayfun(@(x) x.(Property),tracksNA); %(u.^2+v.^2).^(1/2);
if strcmp(PropRange,'auto')
    clear PropRange
    PropRange=[min(intensity,[],"omitnan") max(intensity,[], "omitnan")];
end
vColor=floor(scaleContrast(intensity,PropRange,[1 nColors]));
vColor(vColor<1)=1;
vColor(vColor>nColors)=nColors;
vIndex= unique(vColor);
% Removing NaNs
vIndex = vIndex(~isnan(vIndex));
iiformat = ['%.' '3' 'd'];
hold on
if numel(tracksNA)==1
    s = struct2table(tracksNA,'AsArray',true);
else
    s = struct2table(tracksNA);
end
refineFAID = s.refineFAID;
% Create array of quiverplots
for i=1:numel(vIndex)
    idx = find(vColor==vIndex(i));

    xmat = cell2mat(arrayfun(@(x) x.xCoord(1:iFrame),tracksNA(idx),'UniformOutput',false));
    ymat = cell2mat(arrayfun(@(x) x.yCoord(1:iFrame),tracksNA(idx),'UniformOutput',false));
    if size(xmat,2)==1
        h{i}=plot(xmat',ymat','o','Color',Colormap(vIndex(i),:),'MarkerSize',markerSize);
    else
        h{i}=plot(xmat',ymat','Color',Colormap(vIndex(i),:),'LineWidth',lineWidth);
        plot(xmat(:,end)',ymat(:,end)','o','Color',Colormap(vIndex(i),:),'LineWidth',lineWidth/2,'MarkerSize',markerSize);
    end

    % focal adhesion segmentation
    labelTifPath = [p.OutputDirectory filesep 'labelTifs'];
    maskAdhesion = imread(strcat(labelTifPath,'/label',num2str(iFrame,iiformat),'.tif'));
    labelAdhesion = bwlabel(maskAdhesion,4);
    maxLabel=max(labelAdhesion(:));
    adhBound = cell(maxLabel,1);
    for ii=1:maxLabel
        curAdhBound = bwboundaries(labelAdhesion==ii,4,'noholes');
        adhBound{ii} = curAdhBound{1}; % strongly assumes each has only one boundary
    end

    lengthFAs = cellfun(@length, refineFAID);
    longEnoughFAIDs = lengthFAs >= iFrame;
    idAdhCur = zeros(size(lengthFAs));
    idAdhCur(longEnoughFAIDs) = cellfun(@(x) x(iFrame),refineFAID(longEnoughFAIDs)); 


    %From idx, find relevant refineFAID
    idMappedCur = idAdhCur(idx);
    idMappedCur = idMappedCur(~isnan(idMappedCur) & (idMappedCur>0));
    arrayfun(@(x) plot(x{1}(:,2),x{1}(:,1), ...
        'Color',Colormap(vIndex(i),:), 'LineWidth', 0.5),adhBound(idMappedCur))
end
