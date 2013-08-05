edgeId=4; %5
cellId=6;
plotAll=1;

close all

display('This script works only for edges that exist for all frames!')
if length(groupedNetworks)>1
    display('This script works only when groupedNetworks contains only one network!')
    return;
end

goodCellSet=findCells(groupedNetworks,'kPa',35,'deg',[2 3 4 5 6 7],'errs',0);

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

tmin=1;
tmax=numFrames;

if plotAll
    % plot the time courses:
    figure()
    plot(tmin:tmax,fi(:,1),'b');
    hold on;
    plot(tmin:tmax,-fi_tot(:,1),'r');
    title('fi_x [b] vs. fi_{tot,x} [r]')
    hold off;
    box on;
    set(gca,'fontsize',20,'LineWidth',2)
    
    
    figure()
    plot((tmin:tmax)-tmin+1,fi(tmin:tmax,1),'b');
    hold on;
    plot((tmin:tmax)-tmin+1,-fi_tot(tmin:tmax,1),'r');
    title('fi_x [b] vs. fi_{tot,x} [r]')
    hold off;
    box on;
    set(gca,'fontsize',20,'LineWidth',2)
    
    % plot the time courses:
    figure()
    plot(tmin:tmax,sqrt(sum(fi.^2,2)),'b');
    hold on;
    plot(tmin:tmax,-fi_tot(:,1),'r');
    title('fi_x [b] vs. fi_{tot,x} [r]')
    hold off;
    box on;
    set(gca,'fontsize',20,'LineWidth',2)
    
    % plot the time courses:
    figure()
    plot(tmin:tmax,fi(:,2),'b');
    hold on;
    plot(tmin:tmax,-fi_tot(:,2),'r');
    title('fi_y [b] vs. fi_{tot,y} [r]')
    hold off;
    box on;
    set(gca,'fontsize',20,'LineWidth',2)
    
    
    % plot the time courses:
    figure()
    plot(tmin:tmax,fi(:,1),'b');
    hold on;
    plot(tmin:tmax,-f_res(:,1),'k');
    title('fi_x [b] vs. f_{res,x} [k]')
    hold off;
    box on;
    set(gca,'fontsize',20,'LineWidth',2)
    
    figure()
    plot((tmin:tmax)-tmin+1,fi(tmin:tmax,1),'b');
    hold on;
    plot((tmin:tmax)-tmin+1,-f_res(tmin:tmax,1),'k');
    title('fi_x [b] vs. f_{res,x} [k]')
    hold off;
    box on;
    set(gca,'fontsize',20,'LineWidth',2)
    
    % plot the time courses:
    figure()
    plot(tmin:tmax,fi(:,2),'b');
    hold on;
    plot(tmin:tmax,-f_res(:,2),'k');
    title('fi_y [b] vs. f_{res,y} [k]')
    hold off;
    box on;
    set(gca,'fontsize',20,'LineWidth',2)
    
end


marker=['r','b','m','c','g','y','k'];
label=1;
h=[];
figure()
for iEdge=1:length(horzcat(corrSets.edge.edgeId))
    currh=plot(tmin:tmax, sqrt(sum(corrSets.edge(iEdge).fc.^2,2)), marker(mod(iEdge,7)+1));
    h=horzcat(h,currh(1)); % 1 for marker, 2 for colored line;
    M{label}=['edgeId=',num2str(corrSets.edge(iEdge).edgeId)];
    label=label+1;
    hold on;
end
plot([ 56  56],[0 1000],'--k');
plot([150 150],[0 1000],'--k');
title(['force magnitude of all individual edges connected to cell: ',num2str(cellId)])
xlim([0 numFrames])
ylim([0 110])
legend(h,M);
box on;
set(gca,'fontsize',20,'LineWidth',2)
hold off
clear M



allCellSets=findCells(groupedNetworks,'kPa',35,'deg',[2 3 4 5 6 7],'errs',0);
cellList=[1 2 3 4 5 6];


figure()
label=1;
h=[];
elE_vals_tot=zeros(numFrames,1);
for cellId=cellList
    % take only the cell of interest:
    cellSetId=find([allCellSets.cellId]==cellId);
    goodCellSet=allCellSets(cellSetId);
    [elE_vals]=collectCellValues(groupedNetworks,goodCellSet,'elE');
    elE_vals_tot=nansum([elE_vals_tot elE_vals],2);
    currh=plot(tmin:tmax, elE_vals,marker(mod(cellId,7)+1));
    h=horzcat(h,currh(1)); % 1 for marker, 2 for colored line;
    M{label}=['cellId=',num2str(cellId)];
    label=label+1;
    hold on;
end
plot([ 56  56],[0 1],'--k');
plot([150 150],[0 1],'--k');
title(['elastic energy for cell: ',num2str(cellList)])
xlim([0 numFrames])
ylim([0 0.08])
legend(h,M);
box on;
set(gca,'fontsize',20,'LineWidth',2)
hold off;
clear M

figure()
label=1;
h=[];
sumFi_vals_tot=zeros(numFrames,1);
for cellId=cellList
    % take only the cell of interest:
    cellSetId=find([allCellSets.cellId]==cellId);
    goodCellSet=allCellSets(cellSetId);
    [sumFi_vals]=collectCellValues(groupedNetworks,goodCellSet,'sumFi');
    sumFi_vals_tot=nansum([sumFi_vals_tot sumFi_vals],2);
    currh=plot(tmin:tmax, sumFi_vals,marker(mod(cellId,7)+1));
    h=horzcat(h,currh(1)); % 1 for marker, 2 for colored line;
    M{label}=['cellId=',num2str(cellId)];
    label=label+1;
    hold on;
end
% since counting twice, we have to multiply by 1/2!
sumFi_vals_tot=1/2*sumFi_vals_tot;
plot([ 56  56],[0 1000],'--k');
plot([150 150],[0 1000],'--k');
title(['Sum of all cell-cell forces acting on cell: ',num2str(cellList)])
xlim([0 numFrames])
ylim([0 215])
legend(h,M);
box on;
set(gca,'fontsize',20,'LineWidth',2)
hold off;
clear M

figure()
label=1;
h=[];
resF_vals_tot=zeros(numFrames,1);
for cellId=cellList
    % take only the cell of interest:
    cellSetId=find([allCellSets.cellId]==cellId);
    goodCellSet=allCellSets(cellSetId);
    [resF_vals]=collectCellValues(groupedNetworks,goodCellSet,'resF');
    resF_vals_tot=nansum([resF_vals_tot sqrt(sum(resF_vals.^2,2))],2);
    currh=plot(tmin:tmax, sqrt(sum(resF_vals.^2,2)),marker(mod(cellId,7)+1));
    h=horzcat(h,currh(1)); % 1 for marker, 2 for colored line;
    M{label}=['cellId=',num2str(cellId)];
    label=label+1;
    hold on;
end
plot([ 56  56],[0 1000],'--k');
plot([150 150],[0 1000],'--k');
title(['residual TF magnitude of cell: ',num2str(cellList)])
xlim([0 numFrames])
ylim([0 145])
legend(h,M);
box on;
set(gca,'fontsize',20,'LineWidth',2)
hold off;
clear M


figure()
plot(tmin+1:tmax, elE_vals_tot(tmin+1:tmax));
hold on;
plot([ 56  56],[0 1],'--k');
plot([150 150],[0 1],'--k');
title('Total Elastic Energy of whole cell cluster: ')
xlim([0 numFrames])
ylim([0 0.25])
box on;
set(gca,'fontsize',20,'LineWidth',2)
hold off;


figure()
plot(tmin+1:tmax, sumFi_vals_tot(tmin+1:tmax));
hold on;
plot([ 56  56],[0 1500],'--k');
plot([150 150],[0 1500],'--k');
title('Sum of all cell-cell forces measured in the cluster')
xlim([0 numFrames])
ylim([0 500])
box on;
set(gca,'fontsize',20,'LineWidth',2)
hold off;

figure()
plot(tmin+1:tmax, resF_vals_tot(tmin+1:tmax));
hold on;
plot([ 56  56],[0 1500],'--k');
plot([150 150],[0 1500],'--k');
title('Sum of all rsidual forces')
xlim([0 numFrames])
ylim([0 550])
box on;
set(gca,'fontsize',20,'LineWidth',2)
hold off;

trackedNet_cl1=groupedNetworks.cluster{1}.trackedNet;
label=1;
h=[];
numEdges=length(trackedNet_cl1{end}.edge);
numFrames=length(trackedNet_cl1);
for iEdge=1:numEdges
    for iFrame=1:numFrames
        if length(trackedNet_cl1{iFrame}.edge)>=iEdge && ~isempty(trackedNet_cl1{iFrame}.edge{iEdge})
            edgeID{iEdge}.fcmag(iFrame)=sqrt(sum(trackedNet_cl1{iFrame}.edge{iEdge}.fc.^2,2));
        else
            edgeID{iEdge}.fcmag(iFrame)=NaN;
        end
    end
end

figure()
for iEdge=1:numEdges
    currh=plot(tmin:tmax, edgeID{iEdge}.fcmag, marker(mod(iEdge,7)+1));
    h=horzcat(h,currh(1)); % 1 for marker, 2 for colored line;
    M{label}=['edgeId=',num2str(iEdge)];
    label=label+1;
    hold on;
end
plot([ 56  56],[0 1000],'--k');
plot([150 150],[0 1000],'--k');
title('cell-cell forces at all edges')
xlim([0 numFrames])
ylim([0 110])
legend(h,M);
box on;
set(gca,'fontsize',20,'LineWidth',2)
hold off
clear M