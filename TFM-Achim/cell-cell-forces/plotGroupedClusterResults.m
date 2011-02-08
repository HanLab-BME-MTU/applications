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
% plot elastic energy and residual force over the degree.
%**************************************************************************
goodSet=findCells(groupedClusters,'kPa',10,'myo',0,'myoGlb',-1);
[deg_vals,elE_vals,resF_vals]=collectCellValues(groupedClusters,goodSet,'deg','elE','resF');
figure()
boxplot(elE_vals,deg_vals)
title(['Contraction Energy for cells with connectivity: ',num2str(1:max(deg_vals))])
xlabel('Deg of connectivity')
ylabel('Contraction Energy [pJ]')

figure()
boxplot(sqrt(sum(resF_vals.^2,2)),deg_vals)
title(['Residual force for cells with connectivity: ',num2str(1:max(deg_vals))])
xlabel('Deg of connectivity')
ylabel('Residual force [nN]')

%**************************************************************************
% Compare network and cluster analysis:
%**************************************************************************
goodEdgeSet=findEdges(groupedClusters,'errF',50,'errs',0);
[f1_vals,f2_vals,fc1_vals]=collectEdgeValues(groupedClusters,goodEdgeSet,'f1','f2','fc1');
fnet_mag = sqrt(sum((0.5*(f1_vals-f2_vals)).^2,2));
fc_mag   = sqrt(sum(fc1_vals.^2,2));
figure()
plot(fc_mag,fnet_mag,'*');
xlim([min([fc_mag;fnet_mag])  max([fc_mag;fnet_mag])])
ylim([min([fc_mag;fnet_mag])  max([fc_mag;fnet_mag])])
xlabel('fc [nN]')
ylabel('fnet [nN]')

%**************************************************************************
% plot the interfacial force in depdence of pair degree of connectivity:
%**************************************************************************
goodEdgeSet=findEdges(groupedClusters,'deg',[1 2 3 4 5],'myo',[0 1 -1],'type',{'myoIIA_hp93';'myoIIA_hp94';'myoIIB_hp103'},'errF',50,'errs',0);
[deg_vals,lgth_vals,f1_vals,f2_vals,fc1_vals,fc2_vals]=collectEdgeValues(groupedClusters,goodEdgeSet,'deg','lgth','f1','f2','fc1','fc2');
deg_vals_sorted=sort(deg_vals,2);
fc1_mag = sqrt(sum(fc1_vals.^2,2));
figure()
boxplot(fc1_mag,num2str(deg_vals_sorted))

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
[corrResults]=calCorrResults(corrSets);
