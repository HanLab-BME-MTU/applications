function []=plotGroupedClusterResults(groupedClusters)
if nargin<1 || isempty(groupedClusters)
    try
        load('groupedClusters.mat')
    catch
        display('Couldn''t find any groupedClusters... and stoped!')
        return;
    end    
end

%close all;

onlyCorr=0;
if ~onlyCorr
%**************************************************************************
% Compare network and cluster analysis:
%**************************************************************************
goodEdgeSet=findEdges(groupedClusters);%,'errF',500,'errs',0);
[f1_vals,f2_vals,fc1_vals]=collectEdgeValues(groupedClusters,goodEdgeSet,'f1','f2','fc1');
fnet_mag = sqrt(sum((0.5*(f1_vals-f2_vals)).^2,2));
fc_mag   = sqrt(sum(fc1_vals.^2,2));

figure()
% checked with part7
glbMin=min([fc_mag;fnet_mag]);
glbMax=max([fc_mag;fnet_mag]);
plot(fc_mag,fnet_mag,'*');
hold on
plot([0 glbMax],[0 glbMax],'--k')
xlim([glbMin glbMax])
ylim([glbMin glbMax])
xlabel('fc [nN]')
ylabel('fnet [nN]')
hold off

% Do the angular deviation.

%**************************************************************************
% plot elastic energy and residual force over the degree.
%**************************************************************************
goodCellSet=findCells(groupedClusters,'kPa',8,'myo',1,'myoGlb',[-1 0 1],'errF',Inf,'errs',0);
[deg_vals,elE_vals,sumFmag_vals,resF_vals,sumFi_vals,sumLi_vals]=collectCellValues(groupedClusters,goodCellSet,'deg','elE','sumFmag','resF','sumFi','sumLi');


figure()
% for myosin cells (at least on 8kPa) this ration decreases more quickly
% with the node deg. of connectivity.
boxplot(sqrt(sum(resF_vals.^2,2))./sumFi_vals,deg_vals,'notch','on')
title(['magnitude of residual force / Sum of interfacial forces: ',num2str(1:max(deg_vals))])
xlabel('Deg of connectivity')
ylabel('|residual force| / Sum |interfacial forces|[1]')


figure()
% for myosin cells (at least on 8kPa) this ratio increases more quickly
% with the node deg. of connectivity, this indicates that forces are
% communicted through the cell.
boxplot(sumFi_vals./sumFmag_vals,deg_vals,'notch','on')
title(['Sum of interfacial forces / sum magnitude of traction forces: ',num2str(1:max(deg_vals))])
xlabel('Deg of connectivity')
ylabel(' Sum |interfacial forces| /  sum |traction forces| [1]')


figure()
% checked with part7
boxplot(sumFmag_vals,deg_vals,'notch','on')
title(['Sum of traction force magnitude of cells with connectivity: ',num2str(1:max(deg_vals))])
xlabel('Deg of connectivity')
ylabel('Sum of traction force magnitude [nN]')

figure()
% checked with part7
boxplot(elE_vals,deg_vals,'notch','on')
title(['Contraction Energy for cells with connectivity: ',num2str(1:max(deg_vals))])
xlabel('Deg of connectivity')
ylabel('Contraction Energy [pJ]')

figure()
% checked with part7
boxplot(sqrt(sum(resF_vals.^2,2)),deg_vals,'notch','on')
title(['Residual force for cells with connectivity: ',num2str(1:max(deg_vals))])
xlabel('Deg of connectivity')
ylabel('Residual force [nN]')
set(gca,'fontsize',16)

figure()
% checked with part7
boxplot(sqrt(sum(resF_vals.^2,2))./elE_vals,deg_vals,'notch','on')
title(['Residual force per elastic energy for cells with connectivity: ',num2str(1:max(deg_vals))])
xlabel('Deg of connectivity')
ylabel('Residual force / contraction energy [nN/pJ]')

figure()
% checked with part7
boxplot(sumFi_vals,deg_vals,'notch','on')
title(['Sum of interfacial forces for cells with connectivity: ',num2str(1:max(deg_vals))])
xlabel('Deg of connectivity')
ylabel('Sum of interfacial forces [nN]')
set(gca,'fontsize',16)

figure()
% checked with part7
boxplot(sumLi_vals,deg_vals,'notch','on')
title(['Cumulative length of all interfaces for cells with connectivity: ',num2str(1:max(deg_vals))])
xlabel('Deg of connectivity')
ylabel('Cumulative length of all interfaces [um]')
set(gca,'fontsize',16)

figure()
% checked with part7
boxplot(sumFi_vals./sumLi_vals,deg_vals,'notch','on')
title(['Sum of interfacial forces / Cumulative length of all interfaces for cells with connectivity: ',num2str(1:max(deg_vals))])
xlabel('Deg of connectivity')
ylabel('Sum interf. forces / Cum. intf. length [nN/um]')

figure()
% checked with part7
boxplot(sumFi_vals./elE_vals,deg_vals,'notch','on')
title(['Sum of interfacial forces / el. energy for cells with connectivity: ',num2str(1:max(deg_vals))])
xlabel('Deg of connectivity')
ylabel('Sum of interfacial forces / el. energy [nN/pJ]')
set(gca,'fontsize',16)


figure()
% checked with part7: same but could be colored according to the deg. of
% connectivity as in part7.
plot(sumLi_vals,sumFi_vals,'*')
title(['Correlation: cum. interf. force with cumulative length of all interfaces (any degree): ',num2str(1:max(deg_vals))])
xlabel('Cumulative length of all interfaces [um]')
ylabel('Sum interf. forces [nN]')

figure()
% checked with part7
plot(elE_vals,sqrt(sum(resF_vals.^2,2)),'*')
title(['Correlation: residual force with elastic energy (any degree): ',num2str(1:max(deg_vals))])
xlabel('Elastic energy [pJ]')
ylabel('Residual force [nN]')

figure()
% checked with part7: same but could be colored according to the deg. of
% connectivity as in part7.
plot(elE_vals,sumFi_vals,'*')
title(['Correlation: sum interf. forces with elastic energy (any degree): ',num2str(1:max(deg_vals))])
xlabel('Elastic energy [pJ]')
ylabel('sum interf. forces [nN]')

%**************************************************************************
% Plot sumFi, resF, elE, fi for E=8,35kPa:                                *
%**************************************************************************
% plotResultsForTwoStiff(groupedClusters);


%**************************************************************************
% plot the interfacial force in depdence of pair degree of connectivity:
%**************************************************************************
goodEdgeSet=findEdges(groupedClusters,'kPa',8,'myo',1,'type',{'myoIIA_hp93'},'errF',500,'errs',0);
% goodEdgeSet=findEdges(groupedClusters,'kPa',8,'myo',1,'type',{'myoIIA_hp93';'myoIIA_hp94';'myoIIB_hp103'},'errF',500,'errs',0);
% goodEdgeSet=findEdges(groupedClusters,'kPa',8,'myo',0,'myoGlb',[0],'errF',500,'errs',0);
[deg_vals,lgth_vals,f1_vals,f2_vals,fc1_vals,fc2_vals]=collectEdgeValues(groupedClusters,goodEdgeSet,'deg','lgth','f1','f2','fc1','fc2');
deg_vals_sorted=sort(deg_vals,2);
fc1_mag = sqrt(sum(fc1_vals.^2,2));

figure()
% checked with part7
degij_fcmag=[deg_vals_sorted fc1_mag];
% first sort according to the first column (the small degree value, since
% presorted above), and then sort once more according to the second degree
% value to obtain a nice ordering of the degree pairs:
degij_fcmag_dbl_sorted=sortrows(degij_fcmag,[1 2]);
fc_mag_dbl_sorted = degij_fcmag_dbl_sorted(:,3);
deg_dbl_sorted    = degij_fcmag_dbl_sorted(:,1:2);
boxplot(fc_mag_dbl_sorted,num2str(deg_dbl_sorted),'notch','on')
title(['Interface force for edges with connectivity: ',num2str(1:max(deg_vals))])
xlabel('Deg of connectivity')
ylabel('Interface force [nN]')

figure()
% replot grouped values:
% make the legend:
deg_dbl_sorted_groups=deg_dbl_sorted(:,1);
label=1;
while ~isempty(deg_dbl_sorted_groups)
    degGroup=deg_dbl_sorted_groups(1);
    checkVec=deg_dbl_sorted_groups==degGroup;
    degCount=sum(checkVec);
    deg_dbl_sorted_groups(checkVec)=[];
    M{label}=['deg: ',num2str(degGroup),'-x; [n= ',num2str(degCount),']'];
    label=label+1;
end
boxplot(fc_mag_dbl_sorted,num2str(deg_dbl_sorted(:,1)),'labels',M,'notch','on')
% this should be compared with the following to make sure that the counting is correct:
% h=boxplot(fc_mag_dbl_sorted,num2str(deg_dbl_sorted(:,1)));
title(['Interface force for edges with connectivity: ',num2str(1:max(deg_vals))])
xlabel('Minimal deg. of connectivity')
ylabel('Interface force [nN]')
clear M


figure()
% checked with part7, but the deg values should be sorted in a nicer way!
boxplot(fc1_mag./lgth_vals,num2str(deg_vals_sorted),'notch','on')
title(['Interface stress for edges with connectivity: ',num2str(1:max(deg_vals))])
xlabel('Deg of connectivity')
ylabel('Interface stress [nN/um]')

figure()
% checked with part7
boxplot(lgth_vals,num2str(deg_vals_sorted),'notch','on')
title(['Interface length for edges with connectivity: ',num2str(1:max(deg_vals))])
xlabel('Deg of connectivity')
ylabel('Interface length [um]')

%**************************************************************************
% Plot ctr-ctr, ctr-myo, myo-myo, for degree 1-1 interfaces for E=8,35kPa:*
%**************************************************************************
plotIntForceDeg11(groupedClusters);


%**************************************************************************
% correlate Ecad intensity and interfacial force:
%**************************************************************************
goodEdgeSet=findEdges(groupedClusters,'kPa',8,'errF',500,'errs',0);
[fc1_vals,Itot_vals,Iavg_vals]=collectEdgeValues(groupedClusters,goodEdgeSet,'fc1','Itot','Iavg');
fc1_mag = sqrt(sum(fc1_vals.^2,2));
figure()
plot(fc1_mag,Itot_vals,'*')
title('Correlation Ecad intensity / cell-cell forces')
xlabel('Cell-cell force magnitude [nN]')
ylabel('Integrated Ecad intensity [a.u.]')
figure()
plot(fc1_mag,Iavg_vals,'o')
title('Correlation Ecad intensity / cell-cell forces')
xlabel('Cell-cell force magnitude [nN]')
ylabel('Average Ecad intensity [a.u.]')

end
normVar=1;
maxLag =10;
errF_val_corr=Inf;
%**************************************************************************
% correlate forces for control cells:
%**************************************************************************
goodCellSet=findCells(groupedClusters,'kPa',8,'deg',[2 3 4],'myo',0,'errF',errF_val_corr,'errs',0);
[corrSets]=collectCellValues(groupedClusters,goodCellSet,'corr');
[corrResults]=calCorrResults(corrSets,maxLag,'usefm',normVar);


%**************************************************************************
% correlate forces for myosin cells:
%**************************************************************************
%goodCellSet=findCells(groupedClusters,'kPa',8,'deg',[2 3 4 5 6 7],'myo',1,'type',{'tln1'},'errF',errF_val_corr,'errs',0);
%goodCellSet=findCells(groupedClusters,'kPa',35,'deg',[2 3 4 5 6 7],'myo',1,'type',{'myoIIB_hp103'},'errF',errF_val_corr,'errs',0);
 goodCellSet=findCells(groupedClusters,'kPa',8,'deg',[2 3 4],'myo',1,'type',{'myoIIA_hp93';'myoIIA_hp94'},'errF',errF_val_corr,'errs',0);
if ~isempty(goodCellSet) && ~isempty(goodCellSet(1).cellId)
    [corrSets]=collectCellValues(groupedClusters,goodCellSet,'corr');
    [corrResults]=calCorrResults(corrSets,maxLag,'usefm',normVar);
else
    display('No myosin cells of this type found!')
end
