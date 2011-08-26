function []=plotIntForceDeg11(groupedClusters)

%**************************************************************************
% Myosin experiments:                                                     *
%**************************************************************************

% select only the myosin IIA experiments only:
expList={'shMYOIIA','shMYOIIa'};
[groupedClustersMyoIIA]=preSelectForExp(groupedClusters,expList);

% Plot for degree 1-1 interfaces: myo: 0-0; 1-1; 0-1
goodEdgeSet=findEdges(groupedClustersMyoIIA,'kPa',8,'deg',1,'myo',0,'type',{'ctr';'ctrl';'tln1';'';'myoIIA_hp93';'myoIIA_hp94'},'errs',0);
[lgth_vals,f1_vals,f2_vals,fc1_vals,fc2_vals]=collectEdgeValues(groupedClustersMyoIIA,goodEdgeSet,'lgth','f1','f2','fc1','fc2');
%deg_vals_sorted=sort(deg_vals,2);
fc1_mag_ctr_ctr_4myo = sqrt(sum(fc1_vals.^2,2));
ctr_ctr_median_4myo= median(fc1_mag_ctr_ctr_4myo);
fc1_mag_ctr_ctr_4myo = fc1_mag_ctr_ctr_4myo/ctr_ctr_median_4myo;

fc1_mag_ctr_myo=[];
goodEdgeSet=findEdges(groupedClustersMyoIIA,'kPa',8,'deg',1,'myo',-1,'type',{'myoIIA_hp93';'myoIIA_hp94'},'errs',0);
[lgth_vals,f1_vals,f2_vals,fc1_vals,fc2_vals]=collectEdgeValues(groupedClustersMyoIIA,goodEdgeSet,'lgth','f1','f2','fc1','fc2');
%deg_vals_sorted=sort(deg_vals,2);
fc1_mag_ctr_myo = sqrt(sum(fc1_vals.^2,2))/ctr_ctr_median_4myo;

fc1_mag_myo_myo=[];
goodEdgeSet=findEdges(groupedClustersMyoIIA,'kPa',8,'deg',1,'myo',1,'type',{'myoIIA_hp93';'myoIIA_hp94'},'errs',0);
[lgth_vals,f1_vals,f2_vals,fc1_vals,fc2_vals]=collectEdgeValues(groupedClustersMyoIIA,goodEdgeSet,'lgth','f1','f2','fc1','fc2');
%deg_vals_sorted=sort(deg_vals,2);
fc1_mag_myo_myo = sqrt(sum(fc1_vals.^2,2))/ctr_ctr_median_4myo;

%**************************************************************************
% Talin experiments:                                                      *
%**************************************************************************

expList={'siTLN1'};
[groupedClustersTln1]=preSelectForExp(groupedClusters,expList);

% Plot for degree 1-1 interfaces: myo: 0-0; 1-1; 0-1
%goodEdgeSet=findEdges(groupedClusters,'kPa',35,'deg',1,'myo',0,'type',{'ctr';'ctrl';'myoIIA_hp93';'myoIIA_hp94';'myoIIB_hp103'},'errs',0);
%goodEdgeSet=findEdges(groupedClusters,'kPa',8,'deg',1,'myo',0,'type',{'ctr';'ctrl';'tln1';'';'myoIIA_hp93';'myoIIA_hp94'},'errs',0);
goodEdgeSet=findEdges(groupedClustersTln1,'kPa',8,'deg',1,'myo',0,'type',{'tln1'},'errs',0);
[lgth_vals,f1_vals,f2_vals,fc1_vals,fc2_vals]=collectEdgeValues(groupedClustersTln1,goodEdgeSet,'lgth','f1','f2','fc1','fc2');
%deg_vals_sorted=sort(deg_vals,2);
fc1_mag_ctr_ctr_4tln = sqrt(sum(fc1_vals.^2,2));
ctr_ctr_median_4tln  = median(fc1_mag_ctr_ctr_4tln);
fc1_mag_ctr_ctr_4tln = fc1_mag_ctr_ctr_4tln/ctr_ctr_median_4tln;

% goodEdgeSet=findEdges(groupedClusters,'kPa',35,'deg',1,'myo',-1,'type',{'myoIIA_hp93';'myoIIA_hp94';'myoIIB_hp103'},'errs',0);
goodEdgeSet=findEdges(groupedClustersTln1,'kPa',8,'deg',1,'myo',-1,'type',{'tln1'},'errs',0);
[lgth_vals,f1_vals,f2_vals,fc1_vals,fc2_vals]=collectEdgeValues(groupedClustersTln1,goodEdgeSet,'lgth','f1','f2','fc1','fc2');
%deg_vals_sorted=sort(deg_vals,2);
fc1_mag_ctr_tln = sqrt(sum(fc1_vals.^2,2))/ctr_ctr_median_4tln;

% goodEdgeSet=findEdges(groupedClusters,'kPa',35,'deg',1,'myo',1,'type',{'myoIIA_hp93';'myoIIA_hp94';'myoIIB_hp103'},'errs',0);
goodEdgeSet=findEdges(groupedClustersTln1,'kPa',8,'deg',1,'myo',1,'type',{'tln1'},'errs',0);
[lgth_vals,f1_vals,f2_vals,fc1_vals,fc2_vals]=collectEdgeValues(groupedClustersTln1,goodEdgeSet,'lgth','f1','f2','fc1','fc2');
%deg_vals_sorted=sort(deg_vals,2);
fc1_mag_tln_tln = sqrt(sum(fc1_vals.^2,2))/ctr_ctr_median_4tln;

figure()
maxLp1=max([length(fc1_mag_ctr_ctr_4myo),length(fc1_mag_ctr_myo),length(fc1_mag_myo_myo),length(fc1_mag_ctr_ctr_4tln),length(fc1_mag_ctr_tln),length(fc1_mag_tln_tln)])+1;
fc1_mag_ctr_ctr_4myo( end+1:maxLp1,1)=NaN;
fc1_mag_ctr_myo( end+1:maxLp1,1)=NaN;
fc1_mag_myo_myo( end+1:maxLp1,1)=NaN;
fc1_mag_ctr_ctr_4tln(end+1:maxLp1,1)=NaN;
fc1_mag_ctr_tln(end+1:maxLp1,1)=NaN;
fc1_mag_tln_tln(end+1:maxLp1,1)=NaN;
G1={'ctr-ctr-all'; 'ctr-ctr-tln'; 'ctr-myo'; 'ctr-tln'; 'myo-myo'; 'tln-tln'};
G2={'8kPa'       ; '8kPa'   ; '8kPa'   ; '8kPa'   ; '8kPa'   ;   '8kPa'};
%G2={'8kPa'   ; '35kPa'  ; '8kPa'   ; '35kPa'  ; '8kPa'   ; '35kPa'};
boxplot([fc1_mag_ctr_ctr_4myo fc1_mag_ctr_ctr_4tln fc1_mag_ctr_myo fc1_mag_ctr_tln fc1_mag_myo_myo fc1_mag_tln_tln],{G1 G2},'factorgap',[20 0],'color','bg','notch','on','symbol','+')
title('Interface force for edges with connectivity 1-1')
xlabel('Pair composition')
ylabel('Interface force [nN]')
%ylim([0 130])
set(gca,'fontsize',16)