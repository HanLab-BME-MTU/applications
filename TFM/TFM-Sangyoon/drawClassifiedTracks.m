function [htrackLine, htrackCircles] = drawClassifiedTracks(allDataClass,tracksNA,iFrame,h,showAll)
if nargin <5
    showAll = 0;%false;
end
markerSize = 4;
numGroups = 9;
colors = distinguishable_colors(numGroups,'k');
% switching colors between group 6 and 9
tempColor = colors(6,:);
colors(6,:) = colors(9,:);
colors(9,:) = tempColor;

idCurrent = arrayfun(@(x) x.startingFrameExtra<=iFrame & x.endingFrameExtra>=iFrame,tracksNA);

if ischar(allDataClass{1})
    idGroup{1} = strcmp(allDataClass,'Group1') & idCurrent;
    idGroup{2} = strcmp(allDataClass,'Group2') & idCurrent;
    idGroup{3} = strcmp(allDataClass,'Group3') & idCurrent;
    idGroup{4} = strcmp(allDataClass,'Group4') & idCurrent;
    idGroup{5} = strcmp(allDataClass,'Group5') & idCurrent;
    idGroup{6} = strcmp(allDataClass,'Group6') & idCurrent;
    idGroup{7} = strcmp(allDataClass,'Group7') & idCurrent;
    idGroup{8} = strcmp(allDataClass,'Group8') & idCurrent;
    idGroup{9} = strcmp(allDataClass,'Group9') & idCurrent;
else
    idGroup = allDataClass; % in case this is a cell format
end

if nargin<4 || isempty(h)
    h=gca;
end
switch showAll
    case 10
        iGroups = 1:numGroups;
    case 0
        iGroups = setdiff(1:numGroups,6);
    otherwise
        iGroups = showAll;
end
% if showNoise
%     iGroups = 1:numGroups;
% else
%     iGroups = setdiff(1:numGroups,6);
% end
for ii=iGroups
    if nargin>2
        htrackLine{ii} = arrayfun(@(x) plot(h, x.xCoord(1:iFrame),x.yCoord(1:iFrame),'Color',colors(ii,:)),tracksNA(idGroup{ii}),'UniformOutput',false);
        htrackCircles{ii} = arrayfun(@(x) plot(h, x.xCoord(iFrame),x.yCoord(iFrame),'o','MarkerSize',markerSize,'Color',colors(ii,:)),tracksNA(idGroup{ii}),'UniformOutput',false);
    else
        htrackLine{ii} = arrayfun(@(x) plot(h, x.xCoord,x.yCoord,'Color',colors(ii,:)),tracksNA(idGroup{ii}),'UniformOutput',false);
        htrackCircles{ii} = arrayfun(@(x) plot(h, x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'o','MarkerSize',markerSize,'Color',colors(ii,:)),tracksNA(idGroup{ii}),'UniformOutput',false);
    end    
end
% legend([htrackG1{1} htrackG2{1} htrackG3{1} htrackG4{1} htrackG5{1} htrackG6{1} htrackG7{1} htrackG8{1} htrackG9{1}],{'G1:turn-over','G2:maturing','G3:moving along protruding edge',...
%     'G4:retracting','G5:stable at the edge','G6:noise or very transient','G7:adhesions at stalling edge','G8:strong stable adhesion', 'G9:weak stable adhesion inside'},'TextColor','w','Location','best')
% legend([htrackG1{1} htrackG2{1} htrackG3{1} htrackG4{1} htrackG5{1} htrackG6{1}],{'G1:turn-over','G2:maturing','G3:moving along protruding edge',...
%     'G4:retracting','G5:stable at the edge','G6:noise or very transient'},'TextColor','w','Location','best')
% legend('boxoff')
