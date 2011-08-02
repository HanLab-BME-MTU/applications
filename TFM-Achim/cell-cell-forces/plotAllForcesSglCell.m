edgeId=5; %5
cellId=3;

close all

display('This script works only for edges that exist for all frames!')
if length(groupedNetworks)>1
    display('This script works only when groupedNetworks contains only one network!')
    return;
end

goodCellSet=findCells(groupedNetworks,'kPa',8,'deg',[2 3 4 5 6 7],'myo',1,'type',{'tln1'},'errs',0);

% take only the cell of interest:
cellSetId=find([goodCellSet.cellId]==cellId);
goodCellSet=goodCellSet(cellSetId);
[corrSets]=collectCellValues(groupedNetworks,goodCellSet,'corr');
numFrames=length(corrSets.t);

% take only the edge of interest:
corrSetEdgeId=find([corrSets.edge.edgeId]==edgeId);

fi     = corrSets.edge(corrSetEdgeId).fc;
f_res  = corrSets.resF;
fi_tot = zeros(size(fi));
for edges2sum=find([corrSets.edge.edgeId]~=edgeId)
    fi_tot(:,1) = nansum([fi_tot(:,1) corrSets.edge(edges2sum).fc(:,1)],2);
    fi_tot(:,2) = nansum([fi_tot(:,2) corrSets.edge(edges2sum).fc(:,2)],2);
end
checkVec=sum(fi_tot,2)==0;
fi_tot(checkVec,1)=NaN;
fi_tot(checkVec,2)=NaN;

% plot the time courses:
figure()
plot(1:numFrames,fi(:,1),'b');
hold on;
plot(1:numFrames,-fi_tot(:,1),'r');
title('fi_x [b] vs. fi_{tot,x} [r]')
hold off;
box on;
set(gca,'fontsize',20,'LineWidth',2)

tmin=35;
figure()
plot((tmin:numFrames)-tmin+1,fi(tmin:numFrames,1),'b');
hold on;
plot((tmin:numFrames)-tmin+1,-fi_tot(tmin:numFrames,1),'r');
title('fi_x [b] vs. fi_{tot,x} [r]')
hold off;
box on;
set(gca,'fontsize',20,'LineWidth',2)

RHO = corr(fi(35:177,1),fi_tot(35:177,1))


% plot the time courses:
figure()
plot(1:numFrames,fi(:,2),'b');
hold on;
plot(1:numFrames,-fi_tot(:,2),'r');
title('fi_y [b] vs. fi_{tot,y} [r]')
hold off;
box on;
set(gca,'fontsize',20,'LineWidth',2)


% plot the time courses:
figure()
plot(1:numFrames,fi(:,1),'b');
hold on;
plot(1:numFrames,-f_res(:,1),'k');
title('fi_x [b] vs. f_{res,x} [k]')
hold off;
box on;
set(gca,'fontsize',20,'LineWidth',2)

figure()
plot((tmin:numFrames)-tmin+1,fi(tmin:numFrames,1),'b');
hold on;
plot((tmin:numFrames)-tmin+1,-f_res(tmin:numFrames,1),'k');
title('fi_x [b] vs. f_{res,x} [k]')
hold off;
box on;
set(gca,'fontsize',20,'LineWidth',2)

RHO = corr(fi(35:177,1),f_res(35:177,1))


% plot the time courses:
figure()
plot(1:numFrames,fi(:,2),'b');
hold on;
plot(1:numFrames,-f_res(:,2),'k');
title('fi_y [b] vs. f_{res,y} [k]')
hold off;
box on;
set(gca,'fontsize',20,'LineWidth',2)