function [htrackG] = drawSelectedTracks(tracksNA,idSelected,iFrame,h)
markerSize = 7;
numGroups = 9;
colors = distinguishable_colors(numGroups,'k');
% switching colors between group 6 and 9
tempColor = colors(6,:);
colors(6,:) = colors(9,:);
colors(9,:) = tempColor;

idCurrent = find(arrayfun(@(x) x.startingFrameExtra<=iFrame & x.endingFrameExtra>=iFrame,tracksNA));

idGroup{1} = intersect(idSelected.idGroup1Selected', idCurrent);
idGroup{2} = intersect(idSelected.idGroup2Selected', idCurrent);
idGroup{3} = intersect(idSelected.idGroup3Selected', idCurrent);
idGroup{4} = intersect(idSelected.idGroup4Selected', idCurrent);
idGroup{5} = intersect(idSelected.idGroup5Selected', idCurrent);
idGroup{6} = intersect(idSelected.idGroup6Selected', idCurrent);
idGroup{7} = intersect(idSelected.idGroup7Selected', idCurrent);
idGroup{8} = intersect(idSelected.idGroup8Selected', idCurrent);
idGroup{9} = intersect(idSelected.idGroup9Selected', idCurrent);

if nargin<3 || isempty(h)
    h=gca;
end
% if showNoise
iGroups = 1:numGroups;
% else
%     iGroups = setdiff(1:numGroups,6);
% end
htrackG=cell(numGroups,1);
for ii=iGroups
    xmat = cell2mat(arrayfun(@(x) x.xCoord(1:iFrame),tracksNA(idGroup{ii}),'UniformOutput',false));
    ymat = cell2mat(arrayfun(@(x) x.yCoord(1:iFrame),tracksNA(idGroup{ii}),'UniformOutput',false));
    if size(xmat,2)==1
        htrackG{ii} = plot(xmat',ymat','.','Color',colors(ii,:),'MarkerSize',markerSize);
    else
        htrackG{ii} = plot(xmat',ymat','Color',colors(ii,:));
    end
end
% legend([htrackG{1} htrackG{2} htrackG{3} htrackG{4} htrackG{5} htrackG{6} htrackG{7} htrackG{8} htrackG{9}],{'G1:turn-over','G2:maturing','G3:moving along protruding edge',...
%     'G4:retracting','G5:stable at the edge','G6:noise or very transient','G7:adhesions at stalling edge','G8:strong stable adhesion', 'G9:weak stable adhesion inside'},'TextColor','w','Location','best')
% legend([htrackG{1} htrackG{1} htrackG{1} htrackG{1} htrackG{1} htrackG{1}],{'G1:turn-over','G2:maturing','G3:moving along protruding edge',...
%     'G4:retracting','G5:stable at the edge','G6:noise or very transient'},'TextColor','w','Location','best')
% legend('boxoff')