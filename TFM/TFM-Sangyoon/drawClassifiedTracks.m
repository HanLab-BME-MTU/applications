function [] = drawClassifiedTracks(allDataClass,tracksNA,iFrame)

colors = distinguishable_colors(9,'k');
idGroup1 = strcmp(allDataClass,'Group1');
idGroup2 = strcmp(allDataClass,'Group2');
idGroup3 = strcmp(allDataClass,'Group3');
idGroup4 = strcmp(allDataClass,'Group4');
idGroup5 = strcmp(allDataClass,'Group5');
idGroup6 = strcmp(allDataClass,'Group6');
idGroup7 = strcmp(allDataClass,'Group7');
idGroup8 = strcmp(allDataClass,'Group8');
idGroup9 = strcmp(allDataClass,'Group9');
if nargin>2
    htrackG1=arrayfun(@(x) plot(x.xCoord(1:iFrame),x.yCoord(1:iFrame),'Color',colors(1,:)),tracksNA(idGroup1),'UniformOutput',false);
    htrackG2=arrayfun(@(x) plot(x.xCoord(1:iFrame),x.yCoord(1:iFrame),'Color',colors(2,:)),tracksNA(idGroup2),'UniformOutput',false);
    htrackG3=arrayfun(@(x) plot(x.xCoord(1:iFrame),x.yCoord(1:iFrame),'Color',colors(3,:)),tracksNA(idGroup3),'UniformOutput',false);
    htrackG4=arrayfun(@(x) plot(x.xCoord(1:iFrame),x.yCoord(1:iFrame),'Color',colors(4,:)),tracksNA(idGroup4),'UniformOutput',false);
    htrackG5=arrayfun(@(x) plot(x.xCoord(1:iFrame),x.yCoord(1:iFrame),'Color',colors(5,:)),tracksNA(idGroup5),'UniformOutput',false);
    htrackG6=arrayfun(@(x) plot(x.xCoord(1:iFrame),x.yCoord(1:iFrame),'Color',colors(6,:)),tracksNA(idGroup6),'UniformOutput',false);
    htrackG7=arrayfun(@(x) plot(x.xCoord(1:iFrame),x.yCoord(1:iFrame),'Color',colors(7,:)),tracksNA(idGroup7),'UniformOutput',false);
    htrackG8=arrayfun(@(x) plot(x.xCoord(1:iFrame),x.yCoord(1:iFrame),'Color',colors(8,:)),tracksNA(idGroup8),'UniformOutput',false);
    htrackG9=arrayfun(@(x) plot(x.xCoord(1:iFrame),x.yCoord(1:iFrame),'Color',colors(9,:)),tracksNA(idGroup9),'UniformOutput',false);
    arrayfun(@(x) plot(x.xCoord(iFrame),x.yCoord(iFrame),'o','Color',colors(1,:)),tracksNA(idGroup1));
    arrayfun(@(x) plot(x.xCoord(iFrame),x.yCoord(iFrame),'o','Color',colors(2,:)),tracksNA(idGroup2));
    arrayfun(@(x) plot(x.xCoord(iFrame),x.yCoord(iFrame),'o','Color',colors(3,:)),tracksNA(idGroup3));
    arrayfun(@(x) plot(x.xCoord(iFrame),x.yCoord(iFrame),'o','Color',colors(4,:)),tracksNA(idGroup4));
    arrayfun(@(x) plot(x.xCoord(iFrame),x.yCoord(iFrame),'o','Color',colors(5,:)),tracksNA(idGroup5));
    arrayfun(@(x) plot(x.xCoord(iFrame),x.yCoord(iFrame),'o','Color',colors(6,:)),tracksNA(idGroup6));
    arrayfun(@(x) plot(x.xCoord(iFrame),x.yCoord(iFrame),'o','Color',colors(7,:)),tracksNA(idGroup7));
    arrayfun(@(x) plot(x.xCoord(iFrame),x.yCoord(iFrame),'o','Color',colors(8,:)),tracksNA(idGroup8));
    arrayfun(@(x) plot(x.xCoord(iFrame),x.yCoord(iFrame),'o','Color',colors(9,:)),tracksNA(idGroup9));
else
    htrackG1=arrayfun(@(x) plot(x.xCoord,x.yCoord,'Color',colors(1,:)),tracksNA(idGroup1),'UniformOutput',false);
    htrackG2=arrayfun(@(x) plot(x.xCoord,x.yCoord,'Color',colors(2,:)),tracksNA(idGroup2),'UniformOutput',false);
    htrackG3=arrayfun(@(x) plot(x.xCoord,x.yCoord,'Color',colors(3,:)),tracksNA(idGroup3),'UniformOutput',false);
    htrackG4=arrayfun(@(x) plot(x.xCoord,x.yCoord,'Color',colors(4,:)),tracksNA(idGroup4),'UniformOutput',false);
    htrackG5=arrayfun(@(x) plot(x.xCoord,x.yCoord,'Color',colors(5,:)),tracksNA(idGroup5),'UniformOutput',false);
    htrackG6=arrayfun(@(x) plot(x.xCoord,x.yCoord,'Color',colors(6,:)),tracksNA(idGroup6),'UniformOutput',false);
    htrackG7=arrayfun(@(x) plot(x.xCoord,x.yCoord,'Color',colors(7,:)),tracksNA(idGroup7),'UniformOutput',false);
    htrackG8=arrayfun(@(x) plot(x.xCoord,x.yCoord,'Color',colors(8,:)),tracksNA(idGroup8),'UniformOutput',false);
    htrackG9=arrayfun(@(x) plot(x.xCoord,x.yCoord,'Color',colors(9,:)),tracksNA(idGroup9),'UniformOutput',false);
    arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'o','Color',colors(1,:)),tracksNA(idGroup1));
    arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'o','Color',colors(2,:)),tracksNA(idGroup2));
    arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'o','Color',colors(3,:)),tracksNA(idGroup3));
    arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'o','Color',colors(4,:)),tracksNA(idGroup4));
    arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'o','Color',colors(5,:)),tracksNA(idGroup5));
    arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'o','Color',colors(6,:)),tracksNA(idGroup6));
    arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'o','Color',colors(7,:)),tracksNA(idGroup7));
    arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'o','Color',colors(8,:)),tracksNA(idGroup8));
    arrayfun(@(x) plot(x.xCoord(x.endingFrame),x.yCoord(x.endingFrame),'o','Color',colors(9,:)),tracksNA(idGroup9));
end    
% legend([htrackG1{1} htrackG2{1} htrackG3{1} htrackG4{1} htrackG5{1} htrackG6{1} htrackG7{1} htrackG8{1} htrackG9{1}],{'G1:turn-over','G2:maturing','G3:moving along protruding edge',...
%     'G4:retracting','G5:stable at the edge','G6:noise or very transient','G7:adhesions at stalling edge','G8:strong stable adhesion', 'G9:weak stable adhesion inside'},'TextColor','w','Location','best')
% legend([htrackG1{1} htrackG2{1} htrackG3{1} htrackG4{1} htrackG5{1} htrackG6{1}],{'G1:turn-over','G2:maturing','G3:moving along protruding edge',...
%     'G4:retracting','G5:stable at the edge','G6:noise or very transient'},'TextColor','w','Location','best')
% legend('boxoff')
