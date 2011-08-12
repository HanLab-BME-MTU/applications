function []=plotIntForceDeg11(groupedClusters)

% Plot for degree 1-1 interfaces: myo: 0-0; 1-1; 0-1
goodEdgeSet=findEdges(groupedClusters,'kPa',8,'deg',1,'myo',0,'type',{'myoIIA_hp93';'myoIIA_hp94'},'errs',0);
[lgth_vals,f1_vals,f2_vals,fc1_vals,fc2_vals]=collectEdgeValues(groupedClusters,goodEdgeSet,'lgth','f1','f2','fc1','fc2');
%deg_vals_sorted=sort(deg_vals,2);
fc1_mag_ctr_8 = sqrt(sum(fc1_vals.^2,2));

goodEdgeSet=findEdges(groupedClusters,'kPa',8,'deg',1,'myo',-1,'type',{'myoIIA_hp93';'myoIIA_hp94'},'errs',0);
[lgth_vals,f1_vals,f2_vals,fc1_vals,fc2_vals]=collectEdgeValues(groupedClusters,goodEdgeSet,'lgth','f1','f2','fc1','fc2');
%deg_vals_sorted=sort(deg_vals,2);
fc1_mag_mix_8 = sqrt(sum(fc1_vals.^2,2));

goodEdgeSet=findEdges(groupedClusters,'kPa',8,'deg',1,'myo',1,'type',{'myoIIA_hp93';'myoIIA_hp94'},'errs',0);
[lgth_vals,f1_vals,f2_vals,fc1_vals,fc2_vals]=collectEdgeValues(groupedClusters,goodEdgeSet,'lgth','f1','f2','fc1','fc2');
%deg_vals_sorted=sort(deg_vals,2);
fc1_mag_myo_8 = sqrt(sum(fc1_vals.^2,2));

% % Plot for degree 1-1 interfaces: myo: 0-0; 1-1; 0-1
% goodEdgeSet=findEdges(groupedClusters,'kPa',35,'deg',1,'myo',0,'type',{'myoIIA_hp93';'myoIIA_hp94';'myoIIB_hp103'},'errs',0);
% [lgth_vals,f1_vals,f2_vals,fc1_vals,fc2_vals]=collectEdgeValues(groupedClusters,goodEdgeSet,'lgth','f1','f2','fc1','fc2');
% %deg_vals_sorted=sort(deg_vals,2);
% fc1_mag_ctr_35 = sqrt(sum(fc1_vals.^2,2));
% 
% goodEdgeSet=findEdges(groupedClusters,'kPa',35,'deg',1,'myo',-1,'type',{'myoIIA_hp93';'myoIIA_hp94';'myoIIB_hp103'},'errs',0);
% [lgth_vals,f1_vals,f2_vals,fc1_vals,fc2_vals]=collectEdgeValues(groupedClusters,goodEdgeSet,'lgth','f1','f2','fc1','fc2');
% %deg_vals_sorted=sort(deg_vals,2);
% fc1_mag_mix_35 = sqrt(sum(fc1_vals.^2,2));
% 
% goodEdgeSet=findEdges(groupedClusters,'kPa',35,'deg',1,'myo',1,'type',{'myoIIA_hp93';'myoIIA_hp94';'myoIIB_hp103'},'errs',0);
% [lgth_vals,f1_vals,f2_vals,fc1_vals,fc2_vals]=collectEdgeValues(groupedClusters,goodEdgeSet,'lgth','f1','f2','fc1','fc2');
% %deg_vals_sorted=sort(deg_vals,2);
% fc1_mag_myo_35 = sqrt(sum(fc1_vals.^2,2));

figure()
maxLp1=max([length(fc1_mag_ctr_8),length(fc1_mag_mix_8),length(fc1_mag_myo_8),length(fc1_mag_ctr_35),length(fc1_mag_mix_35),length(fc1_mag_myo_35)])+1;
fc1_mag_ctr_8( end+1:maxLp1,1)=NaN;
fc1_mag_mix_8( end+1:maxLp1,1)=NaN;
fc1_mag_myo_8( end+1:maxLp1,1)=NaN;
fc1_mag_ctr_35(end+1:maxLp1,1)=NaN;
fc1_mag_mix_35(end+1:maxLp1,1)=NaN;
fc1_mag_myo_35(end+1:maxLp1,1)=NaN;
G1={'ctr-ctr'; 'ctr-ctr'; 'ctr-myo'; 'ctr-myo'; 'myo-myo'; 'myo-myo'};
G2={'8kPa'   ; '35kPa'  ; '8kPa'   ; '35kPa'  ; '8kPa'   ; '35kPa'};
boxplot([fc1_mag_ctr_8 fc1_mag_ctr_35 fc1_mag_mix_8 fc1_mag_mix_35 fc1_mag_myo_8 fc1_mag_myo_35],{G1 G2},'factorgap',[20 0],'color','bg','notch','on','symbol','+')
title('Interface force for edges with connectivity 1-1')
xlabel('Pair composition')
ylabel('Interface force [nN]')
ylim([0 120])
set(gca,'fontsize',16)