function []=plotResultsForTwoStiff(groupedClusters,degcutoff)
if nargin<2
    degcutoff=5;
else
    degcutoff=Inf;
end
display(['No degrees >=',num2str(degcutoff),' are considered'])

goodCellSet_35=findCells(groupedClusters,'kPa',35,'myo',0,'myoGlb',[0],'errF',Inf,'errs',0);
[deg_vals_35,elE_vals_35,sumFmag_vals_35,resF_vals_35,sumFi_vals_35,sumLi_vals_35]=collectCellValues(groupedClusters,goodCellSet_35,'deg','elE','sumFmag','resF','sumFi','sumLi');
cv=(deg_vals_35>=degcutoff | isnan(deg_vals_35));
deg_vals_35(cv)=[]; elE_vals_35(cv)=[]; sumFmag_vals_35(cv)=[]; resF_vals_35(cv,:)=[]; sumFi_vals_35(cv)=[]; sumLi_vals_35(cv)=[];

goodCellSet_8=findCells(groupedClusters,'kPa',8,'myo',0,'myoGlb',0,'errF',Inf,'errs',0);
[deg_vals_8, elE_vals_8, sumFmag_vals_8, resF_vals_8, sumFi_vals_8, sumLi_vals_8]=collectCellValues(groupedClusters,goodCellSet_8,'deg','elE','sumFmag','resF','sumFi','sumLi');
cv=(deg_vals_8>=degcutoff | isnan(deg_vals_8));
deg_vals_8(cv)=[]; elE_vals_8(cv)=[]; sumFmag_vals_8(cv)=[]; resF_vals_8(cv,:)=[]; sumFi_vals_8(cv)=[]; sumLi_vals_8(cv)=[];

%kill all data with deg >=5


n_8 =hist(deg_vals_8 ,[1:15]);
n_35=hist(deg_vals_35,[1:15]);

for k=1:15
    M_8{k} =['8kPa: ',num2str(n_8(k)),';'];
    M_35{k}=['35kPa: ',num2str(n_35(k))];    
end

maxDeg=(max([deg_vals_35;deg_vals_8]));
G8kPa    =repmat({'8kPa'} ,length(deg_vals_8) ,1);
G35kPa   =repmat({'35kPa'},length(deg_vals_35),1);
% Gfake =repmat({''      ,maxDeg,1);
Gcol=vertcat(G8kPa,G35kPa);

G1=vertcat(deg_vals_8,deg_vals_35);
G2=[M_8(deg_vals_8),M_35(deg_vals_35)]';

figure()
% checked with part7
boxplot([sumLi_vals_8;sumLi_vals_35],{G1 G2},'factorgap',[20 0],'colorgroup',Gcol,'color','br','notch','on','symbol','.')
title(['Cumulative length of all interfaces for cells with connectivity: ',num2str(1:maxDeg)])
xlabel('Deg of connectivity')
ylabel('Cumulative length of all interfaces [um]')
set(gca,'fontsize',20)
box on
ylim([0,115])

figure()
% checked with part7
boxplot([sumFi_vals_8;sumFi_vals_35],{G1 G2},'factorgap',[20 0],'colorgroup',Gcol,'color','br','notch','on','symbol','.')
title(['Sum of interface forces for degrees: ',num2str(1:maxDeg)])
xlabel('Deg of connectivity')
ylabel('sum interface forces [nN]')
set(gca,'fontsize',20)
box on
% ylim([0,275])
ylim([0,355])

figure()
% checked with part7
boxplot(sqrt(sum([resF_vals_8; resF_vals_35].^2,2)),{G1 G2},'factorgap',[20 0],'colorgroup',Gcol,'color','br','notch','on','symbol','.')
title(['Residual force for degrees: ',num2str(1:maxDeg)])
xlabel('Deg of connectivity')
ylabel('Residual force [nN]')
set(gca,'fontsize',20)
box on
%ylim([0,160])
ylim([0,305])

figure()
% checked with part7
boxplot([elE_vals_8; elE_vals_35],{G1 G2},'factorgap',[20 0],'colorgroup',Gcol,'color','br','notch','on','symbol','.')
title(['Contraction Energy for cells with connectivity: ',num2str(1:maxDeg)])
xlabel('Deg of connectivity')
ylabel('Contraction Energy [pJ]')
set(gca,'fontsize',20)
box on
ylim([0,0.3])


figure()
% checked with part7
boxplot([sumFmag_vals_8; sumFmag_vals_35],{G1 G2},'factorgap',[20 0],'colorgroup',Gcol,'color','br','notch','on','symbol','.')
title(['Sum of traction force magnitudes for cells with connectivity: ',num2str(1:maxDeg)])
xlabel('Deg of connectivity')
ylabel('Sum of traction force magnitudes [nN]')
set(gca,'fontsize',20)
ylim([0,1400])
box on

figure()
% checked with part7
boxplot([sumFi_vals_8; sumFi_vals_35]./[sumFmag_vals_8; sumFmag_vals_35],{G1 G2},'factorgap',[20 0],'colorgroup',Gcol,'color','br','notch','on','symbol','.')
title(['ratio: sum Fi /sum |Ft| for cells with connectivity: ',num2str(1:maxDeg)])
xlabel('Deg of connectivity')
ylabel('sum Fi / ratio sum |Ft|')
set(gca,'fontsize',20)
box on
ylim([0,1.5])

figure()
% checked with part7
boxplot([sumFmag_vals_8; sumFmag_vals_35]+[sumFi_vals_8; sumFi_vals_35],{G1 G2},'factorgap',[20 0],'colorgroup',Gcol,'color','br','notch','on','symbol','.')
title(['total contractility for cells with connectivity: ',num2str(1:maxDeg)])
xlabel('Deg of connectivity')
ylabel('total contractility [nN]')
set(gca,'fontsize',20)
ylim([0,1600])
box on

%plot bla
clear M_8 M_35;

goodEdgeSet_35=findEdges(groupedClusters,'kPa',35,'myo',0,'myoGlb',[0],'errF',Inf,'errs',0);
[deg_vals_35,lgth_vals_35,fc1_vals_35,nVec_vals_35]=collectEdgeValues(groupedClusters,goodEdgeSet_35,'deg','lgth','fc1','nVec');
cv=(sum(deg_vals_35>=degcutoff,2)==2);
deg_vals_35(cv,:)=[]; lgth_vals_35(cv)=[]; fc1_vals_35(cv,:)=[]; nVec_vals_35(cv,:)=[];
deg_vals_sorted_35=sort(deg_vals_35,2);
fc1_mag_35        = sqrt(sum(fc1_vals_35.^2,2));
intfStress_35     = fc1_mag_35./lgth_vals_35;
alpha_fc1_nVec_35 = acosd(dot(fc1_vals_35,nVec_vals_35,2)./(fc1_mag_35));
cv=alpha_fc1_nVec_35>90;
alpha_fc1_nVec_35(cv)=180-alpha_fc1_nVec_35(cv);

goodEdgeSet_8=findEdges(groupedClusters,'kPa',8,'myo',0,'myoGlb',0,'errF',Inf,'errs',0);
[deg_vals_8,lgth_vals_8,fc1_vals_8,nVec_vals_8]=collectEdgeValues(groupedClusters,goodEdgeSet_8,'deg','lgth','fc1','nVec');
cv=(sum(deg_vals_8>=degcutoff,2)==2);
deg_vals_8(cv,:)=[]; lgth_vals_8(cv)=[]; fc1_vals_8(cv,:)=[]; nVec_vals_8(cv,:)=[];
deg_vals_sorted_8 =sort(deg_vals_8,2);
fc1_mag_8         = sqrt(sum(fc1_vals_8.^2,2));
intfStress_8      = fc1_mag_8./lgth_vals_8;
alpha_fc1_nVec_8 = acosd(dot(fc1_vals_8,nVec_vals_8,2)./(fc1_mag_8));
%quick cheat:
cv=alpha_fc1_nVec_8>90;
alpha_fc1_nVec_8(cv)=180-alpha_fc1_nVec_8(cv);

% figure()
[deg_dbl_sorted_35,idx_35]  =sortrows(deg_vals_sorted_35,[1 2]);
fc_mag_dbl_sorted_35        =fc1_mag_35(idx_35);
intfStress_dbl_sorted_35    =intfStress_35(idx_35);
alpha_fc1_nVec_dbl_sorted_35=alpha_fc1_nVec_35(idx_35);
% boxplot(fc_mag_dbl_sorted_35,num2str(deg_dbl_sorted_35),'notch','on')
% title(['Interface force for edges with connectivity: ',num2str(1:max(deg_vals_35))])
% xlabel('Deg of connectivity')
% ylabel('Interface force [nN]')

% figure()
deg_dbl_sorted_groups_35=deg_dbl_sorted_35(:,1);
label=1;
while ~isempty(deg_dbl_sorted_groups_35)
    degGroup_35=deg_dbl_sorted_groups_35(1);
    checkVec_35=deg_dbl_sorted_groups_35==degGroup_35;
    degCount_35=sum(checkVec_35);
    deg_dbl_sorted_groups_35(checkVec_35)=[];
    %M_35{label}=['deg: ',num2str(degGroup_35),'-x; [n= ',num2str(degCount_35),']'];
    M_35{label}=['35kPa: ',num2str(degCount_35)];
    label=label+1;
end
% boxplot(fc_mag_dbl_sorted_35,num2str(deg_dbl_sorted_35(:,1)),'labels',M_35,'notch','on')
% title(['Interface force for edges with connectivity: ',num2str(1:max(deg_vals_35))])
% xlabel('Minimal deg. of connectivity')
% ylabel('Interface force [nN]')
%% clear M_35


% figure()
[deg_dbl_sorted_8,idx_8]   =sortrows(deg_vals_sorted_8,[1 2]);
fc_mag_dbl_sorted_8        =fc1_mag_8(idx_8);
intfStress_dbl_sorted_8    =intfStress_8(idx_8);
alpha_fc1_nVec_dbl_sorted_8=alpha_fc1_nVec_8(idx_8);
% boxplot(fc_mag_dbl_sorted_8,num2str(deg_dbl_sorted_8),'notch','on')
% title(['Interface force for edges with connectivity: ',num2str(1:max(deg_vals_8))])
% xlabel('Deg of connectivity')
% ylabel('Interface force [nN]')

% figure()
deg_dbl_sorted_groups_8=deg_dbl_sorted_8(:,1);
label=1;
while ~isempty(deg_dbl_sorted_groups_8)
    degGroup_8=deg_dbl_sorted_groups_8(1);
    checkVec_8=deg_dbl_sorted_groups_8==degGroup_8;
    degCount_8=sum(checkVec_8);
    deg_dbl_sorted_groups_8(checkVec_8)=[];
    M_8{label}=['8kPa: ',num2str(degCount_8),';'];
    label=label+1;
end
% boxplot(fc_mag_dbl_sorted_8,num2str(deg_dbl_sorted_8(:,1)),'labels',M_8,'notch','on')
% title(['Interface force for edges with connectivity: ',num2str(1:max(deg_vals_8))])
% xlabel('Minimal deg. of connectivity')
% ylabel('Interface force [nN]')
% %clear M_8


maxDeg=(max([deg_dbl_sorted_8(:,1);deg_dbl_sorted_35(:,1)]));
G8kPa    =repmat({'8kPa'} ,length(deg_dbl_sorted_8(:,1)) ,1);
G35kPa   =repmat({'35kPa'},length(deg_dbl_sorted_35(:,1)),1);

G1=[deg_dbl_sorted_8(:,1);deg_dbl_sorted_35(:,1)];
Gcol=vertcat(G8kPa,G35kPa);

G2=[M_8(deg_dbl_sorted_8(:,1)),M_35(deg_dbl_sorted_35(:,1))]';

figure()
% plot both together:
% 'labels',{[M_8 M_35]}
boxplot([fc_mag_dbl_sorted_8;fc_mag_dbl_sorted_35],{G1 G2},'factorgap',[20 0],'colorgroup',Gcol,'color','br','notch','on','symbol','.')
% this should be compared with the following to make sure that the counting is correct:
% h=boxplot(fc_mag_dbl_sorted,num2str(deg_dbl_sorted(:,1)));
title(['Interface force for edges with connectivity: ',num2str(1:maxDeg)])
xlabel('Minimal deg. of connectivity')
ylabel('Interface force [nN]')
set(gca,'fontsize',20)
%ylim([0,115])
ylim([0,220])
box on

figure()
% plot both together:
boxplot([intfStress_dbl_sorted_8;intfStress_dbl_sorted_35],{G1 G2},'factorgap',[20 0],'colorgroup',Gcol,'color','br','notch','on','symbol','.')
% this should be compared with the following to make sure that the counting is correct:
% h=boxplot(fc_mag_dbl_sorted,num2str(deg_dbl_sorted(:,1)));
title(['Interface stress for edges with connectivity: ',num2str(1:maxDeg)])
xlabel('Minimal deg. of connectivity')
ylabel('Interface stress [nN/um]')
set(gca,'fontsize',20)
%ylim([0,8])
ylim([0,9])
box on


% % plot the stress:
% figure()
% % checked with part7, but the deg values should be sorted in a nicer way!
% boxplot(intfStress_35,num2str(deg_vals_sorted_35),'notch','on')
% title(['Interface stress for edges with connectivity: ',num2str(1:max(deg_vals_35))])
% xlabel('Deg of connectivity')
% ylabel('Interface stress [nN/um]')
% 
% figure()
% % checked with part7
% boxplot(lgth_vals_35,num2str(deg_vals_sorted_35),'notch','on')
% title(['Interface length for edges with connectivity: ',num2str(1:max(deg_vals_35))])
% xlabel('Deg of connectivity')
% ylabel('Interface length [um]')


clear M_8 M_35;

goodEdgeSet_35=findEdges(groupedClusters,'kPa',35,'myo',0,'myoGlb',[0],'errF',Inf,'errs',0);
[deg_vals_35,lgth_vals_35,fc1_vals_35,nVec_vals_35]=collectEdgeValues(groupedClusters,goodEdgeSet_35,'deg','lgth','fc1','nVec');
cv=(sum(deg_vals_35>=degcutoff,2)==2);
deg_vals_35(cv,:)=[]; lgth_vals_35(cv)=[]; fc1_vals_35(cv,:)=[]; nVec_vals_35(cv,:)=[];
deg_vals_sorted_35=sort(deg_vals_35,2);
fc1_mag_35        = sqrt(sum(fc1_vals_35.^2,2));
alpha_fc1_nVec_35 = acosd(dot(fc1_vals_35,nVec_vals_35,2)./(fc1_mag_35));
cv=alpha_fc1_nVec_35>90;
alpha_fc1_nVec_35(cv)=180-alpha_fc1_nVec_35(cv);

goodEdgeSet_8=findEdges(groupedClusters,'kPa',8,'myo',0,'myoGlb',0,'errF',Inf,'errs',0);
% goodEdgeSet_8=findEdges(groupedClusters,'kPa',8,'errF',30,'errs',0);
[deg_vals_8,lgth_vals_8,fc1_vals_8,nVec_vals_8]=collectEdgeValues(groupedClusters,goodEdgeSet_8,'deg','lgth','fc1','nVec');
cv=(sum(deg_vals_8>=degcutoff,2)==2);
deg_vals_8(cv,:)=[]; lgth_vals_8(cv)=[]; fc1_vals_8(cv,:)=[]; nVec_vals_8(cv,:)=[];
deg_vals_sorted_8 =sort(deg_vals_8,2);
fc1_mag_8         = sqrt(sum(fc1_vals_8.^2,2));
alpha_fc1_nVec_8 = acosd(dot(fc1_vals_8,nVec_vals_8,2)./(fc1_mag_8));
%quick cheat:
cv=alpha_fc1_nVec_8>90;
alpha_fc1_nVec_8(cv)=180-alpha_fc1_nVec_8(cv);


[deg_dbl_sorted_35,idx_35]  =sortrows(deg_vals_sorted_35,[1 2]);
fc_mag_dbl_sorted_35        =fc1_mag_35(idx_35);
alpha_fc1_nVec_dbl_sorted_35=alpha_fc1_nVec_35(idx_35);

deg_dbl_sorted_groups_35=deg_dbl_sorted_35(:,1);
label=1;
while ~isempty(deg_dbl_sorted_groups_35)
    degGroup_35=deg_dbl_sorted_groups_35(1);
    checkVec_35=deg_dbl_sorted_groups_35==degGroup_35;
    degCount_35=sum(checkVec_35);
    deg_dbl_sorted_groups_35(checkVec_35)=[];
    %M_35{label}=['deg: ',num2str(degGroup_35),'-x; [n= ',num2str(degCount_35),']'];
    M_35{label}=['35kPa: ',num2str(degCount_35)];
    label=label+1;
end


[deg_dbl_sorted_8,idx_8]  =sortrows(deg_vals_sorted_8,[1 2]);
fc_mag_dbl_sorted_8        =fc1_mag_8(idx_8);
alpha_fc1_nVec_dbl_sorted_8=alpha_fc1_nVec_8(idx_8);

deg_dbl_sorted_groups_8=deg_dbl_sorted_8(:,1);
label=1;
while ~isempty(deg_dbl_sorted_groups_8)
    degGroup_8=deg_dbl_sorted_groups_8(1);
    checkVec_8=deg_dbl_sorted_groups_8==degGroup_8;
    degCount_8=sum(checkVec_8);
    deg_dbl_sorted_groups_8(checkVec_8)=[];
    %M_8{label}=['deg: ',num2str(degGroup_8),'-x; [n= ',num2str(degCount_8),']'];
    M_8{label}=['8kPa: ',num2str(degCount_8)];
    label=label+1;
end

maxDeg=(max([deg_dbl_sorted_8(:,1);deg_dbl_sorted_35(:,1)]));
G8kPa    =repmat({'8kPa'} ,length(deg_dbl_sorted_8(:,1)) ,1);
G35kPa   =repmat({'35kPa'},length(deg_dbl_sorted_35(:,1)),1);

G1=[deg_dbl_sorted_8(:,1);deg_dbl_sorted_35(:,1)];
Gcol=vertcat(G8kPa,G35kPa);

G2=[M_8(deg_dbl_sorted_8(:,1)),M_35(deg_dbl_sorted_35(:,1))]';

figure()
% plot both together:
boxplot([alpha_fc1_nVec_dbl_sorted_8;alpha_fc1_nVec_dbl_sorted_35],{G1 G2},'factorgap',[20 0],'colorgroup',Gcol,'color','br','notch','on','symbol','.')
% this should be compared with the following to make sure that the counting is correct:
% h=boxplot(fc_mag_dbl_sorted,num2str(deg_dbl_sorted(:,1)));
title(['angle: ',num2str(1:maxDeg)])
xlabel('Minimal deg. of connectivity')
ylabel('angle [deg]')
set(gca,'fontsize',20)
box on