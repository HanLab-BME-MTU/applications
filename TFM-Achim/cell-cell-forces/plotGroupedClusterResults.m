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

%**************************************************************************
% Compare network and cluster analysis:
%**************************************************************************
goodEdgeSet=findEdges(groupedClusters);%,'errF',50,'errs',0);
[f1_vals,f2_vals,fc1_vals]=collectEdgeValues(groupedClusters,goodEdgeSet,'f1','f2','fc1');
fnet_mag = sqrt(sum((0.5*(f1_vals-f2_vals)).^2,2));
fc_mag   = sqrt(sum(fc1_vals.^2,2));

figure()
% checked with part7
plot(fc_mag,fnet_mag,'*');
xlim([min([fc_mag;fnet_mag])  max([fc_mag;fnet_mag])])
ylim([min([fc_mag;fnet_mag])  max([fc_mag;fnet_mag])])
xlabel('fc [nN]')
ylabel('fnet [nN]')

% Do the angular deviation.

%**************************************************************************
% plot elastic energy and residual force over the degree.
%**************************************************************************
goodCellSet=findCells(groupedClusters,'kPa',10,'myo',[0 1],'myoGlb',[-1 0 1]);
[deg_vals,elE_vals,resF_vals,sumFi_vals,sumLi_vals]=collectCellValues(groupedClusters,goodCellSet,'deg','elE','resF','sumFi','sumLi');

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

figure()
% checked with part7
boxplot(sumLi_vals,deg_vals,'notch','on')
title(['Cumulative length of all interfaces for cells with connectivity: ',num2str(1:max(deg_vals))])
xlabel('Deg of connectivity')
ylabel('Cumulative length of all interfaces [um]')

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
ylabel('Residual force [nN]')


%**************************************************************************
% plot the interfacial force in depdence of pair degree of connectivity:
%**************************************************************************
goodEdgeSet=findEdges(groupedClusters,'deg',[1 2 3 4 5],'myo',[0 1 -1],'type',{'myoIIA_hp93';'myoIIA_hp94';'myoIIB_hp103'},'errF',50,'errs',0);
[deg_vals,lgth_vals,f1_vals,f2_vals,fc1_vals,fc2_vals]=collectEdgeValues(groupedClusters,goodEdgeSet,'deg','lgth','f1','f2','fc1','fc2');
deg_vals_sorted=sort(deg_vals,2);
fc1_mag = sqrt(sum(fc1_vals.^2,2));

figure()
% checked with part7
boxplot(fc1_mag,num2str(deg_vals_sorted))
title(['Interface force for edges with connectivity: ',num2str(1:max(deg_vals))])
xlabel('Deg of connectivity')
ylabel('Interface force [nN]')

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
% correlate Ecad intensity and interfacial force:
%**************************************************************************
goodEdgeSet=findEdges(groupedClusters,'errF',50,'errs',0);
[fc1_vals,Itot_vals,Iavg_vals]=collectEdgeValues(groupedClusters,goodEdgeSet,'fc1','Itot','Iavg');
fc1_mag = sqrt(sum(fc1_vals.^2,2));
figure()
plot(fc1_mag,Itot_vals,'*')
figure()
plot(fc1_mag,Iavg_vals,'o')

%**************************************************************************
% correlate forces for control cells:
%**************************************************************************
goodCellSet=findCells(groupedClusters,'deg',[2 3 4 5 6],'myo',0,'errs',0);
[corrSets]=collectCellValues(groupedClusters,goodCellSet,'corr');
[corrResults]=calCorrResults(corrSets);


%**************************************************************************
% correlate forces for myosin cells:
%**************************************************************************
goodCellSet=findCells(groupedClusters,'deg',[2 3 4 5 6],'myo',1,'type',{'myoIIA_hp93';'myoIIA_hp94';'myoIIB_hp103'},'errs',0);
[corrSets]=collectCellValues(groupedClusters,goodCellSet,'corr');
maxLag=1;
[corrResults]=calCorrResults(corrSets,maxLag,'usefm');
