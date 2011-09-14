function []=plotResultsForTwoStiff(groupedClusters)

goodCellSet_35=findCells(groupedClusters,'kPa',35,'myo',0,'myoGlb',[0],'errF',Inf,'errs',0);
[deg_vals_35,elE_vals_35,sumFmag_vals_35,resF_vals_35,sumFi_vals_35,sumLi_vals_35]=collectCellValues(groupedClusters,goodCellSet_35,'deg','elE','sumFmag','resF','sumFi','sumLi');
cv=isnan(deg_vals_35);
deg_vals_35(cv)=[]; elE_vals_35(cv)=[]; sumFmag_vals_35(cv)=[]; resF_vals_35(cv,:)=[]; sumFi_vals_35(cv)=[]; sumLi_vals_35(cv)=[];

goodCellSet_8=findCells(groupedClusters,'kPa',8,'myo',0,'myoGlb',0,'errF',Inf,'errs',0);
[deg_vals_8, elE_vals_8, sumFmag_vals_8, resF_vals_8, sumFi_vals_8, sumLi_vals_8]=collectCellValues(groupedClusters,goodCellSet_8,'deg','elE','sumFmag','resF','sumFi','sumLi');
cv=isnan(deg_vals_8);
deg_vals_8(cv)=[]; elE_vals_8(cv)=[]; sumFmag_vals_8(cv)=[]; resF_vals_8(cv,:)=[]; sumFi_vals_8(cv)=[]; sumLi_vals_8(cv)=[];

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
boxplot([sumFi_vals_8;sumFi_vals_35],{G1 G2},'factorgap',[20 0],'colorgroup',Gcol,'color','bg','notch','on','symbol','+')
title(['Sum of interface forces for degrees: ',num2str(1:maxDeg)])
xlabel('Deg of connectivity')
ylabel('sum interface forces [nN]')
set(gca,'fontsize',20)
box on
ylim([0,275])

figure()
% checked with part7
boxplot(sqrt(sum([resF_vals_8; resF_vals_35].^2,2)),{G1 G2},'factorgap',[20 0],'colorgroup',Gcol,'color','bg','notch','on','symbol','+')
title(['Residual force for degrees: ',num2str(1:maxDeg)])
xlabel('Deg of connectivity')
ylabel('Residual force [nN]')
set(gca,'fontsize',20)
box on
ylim([0,160])

figure()
% checked with part7
boxplot([elE_vals_8; elE_vals_35],{G1 G2},'factorgap',[20 0],'colorgroup',Gcol,'color','bg','notch','on','symbol','+')
title(['Contraction Energy for cells with connectivity: ',num2str(1:maxDeg)])
xlabel('Deg of connectivity')
ylabel('Contraction Energy [pJ]')
set(gca,'fontsize',20)
box on
ylim([0,0.11])

clear M_8 M_35;

goodEdgeSet_35=findEdges(groupedClusters,'kPa',35,'myo',0,'myoGlb',[0],'errF',Inf,'errs',0);
[deg_vals_35,lgth_vals_35,fc1_vals_35,nVec_vals_35]=collectEdgeValues(groupedClusters,goodEdgeSet_35,'deg','lgth','fc1','nVec');
deg_vals_sorted_35=sort(deg_vals_35,2);
fc1_mag_35        = sqrt(sum(fc1_vals_35.^2,2));
intfStress_35     = fc1_mag_35./lgth_vals_35;
alpha_fc1_nVec_35 = acosd(dot(fc1_vals_35,nVec_vals_35,2)./(fc1_mag_35));
cv=alpha_fc1_nVec_35>90;
alpha_fc1_nVec_35(cv)=180-alpha_fc1_nVec_35(cv);

goodEdgeSet_8=findEdges(groupedClusters,'kPa',8,'myo',0,'myoGlb',0,'errF',Inf,'errs',0);
[deg_vals_8,lgth_vals_8,fc1_vals_8,nVec_vals_8]=collectEdgeValues(groupedClusters,goodEdgeSet_8,'deg','lgth','fc1','nVec');
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
boxplot(fc_mag_dbl_sorted_8,num2str(deg_dbl_sorted_8(:,1)),'labels',M_8,'notch','on')
title(['Interface force for edges with connectivity: ',num2str(1:max(deg_vals_8))])
xlabel('Minimal deg. of connectivity')
ylabel('Interface force [nN]')
%clear M_8


maxDeg=(max([deg_dbl_sorted_8(:,1);deg_dbl_sorted_35(:,1)]));
G8kPa    =repmat({'8kPa'} ,length(deg_dbl_sorted_8(:,1)) ,1);
G35kPa   =repmat({'35kPa'},length(deg_dbl_sorted_35(:,1)),1);

G1=[deg_dbl_sorted_8(:,1);deg_dbl_sorted_35(:,1)];
Gcol=vertcat(G8kPa,G35kPa);

G2=[M_8(deg_dbl_sorted_8(:,1)),M_35(deg_dbl_sorted_35(:,1))]';

figure()
% plot both together:
% 'labels',{[M_8 M_35]}
boxplot([fc_mag_dbl_sorted_8;fc_mag_dbl_sorted_35],{G1 G2},'factorgap',[20 0],'colorgroup',Gcol,'color','bg','notch','on','symbol','+')
% this should be compared with the following to make sure that the counting is correct:
% h=boxplot(fc_mag_dbl_sorted,num2str(deg_dbl_sorted(:,1)));
title(['Interface force for edges with connectivity: ',num2str(1:maxDeg)])
xlabel('Minimal deg. of connectivity')
ylabel('Interface force [nN]')
set(gca,'fontsize',20)
ylim([0,115])
box on

figure()
% plot both together:
boxplot([intfStress_dbl_sorted_8;intfStress_dbl_sorted_35],{G1 G2},'factorgap',[20 0],'colorgroup',Gcol,'color','bg','notch','on','symbol','+')
% this should be compared with the following to make sure that the counting is correct:
% h=boxplot(fc_mag_dbl_sorted,num2str(deg_dbl_sorted(:,1)));
title(['Interface stress for edges with connectivity: ',num2str(1:maxDeg)])
xlabel('Minimal deg. of connectivity')
ylabel('Interface stress [nN/um]')
set(gca,'fontsize',20)
ylim([0,8])
box on

figure()
% plot both together:
boxplot([alpha_fc1_nVec_dbl_sorted_8;alpha_fc1_nVec_dbl_sorted_35],{G1 G2},'factorgap',[20 0],'colorgroup',Gcol,'color','bg','notch','on','symbol','+')
% this should be compared with the following to make sure that the counting is correct:
% h=boxplot(fc_mag_dbl_sorted,num2str(deg_dbl_sorted(:,1)));
title(['Interface stress for edges with connectivity: ',num2str(1:maxDeg)])
xlabel('Minimal deg. of connectivity')
ylabel('Interface stress [nN/um]')
set(gca,'fontsize',20)
box on


figure()
% checked with part7
degij_alpha_fc_nVec=[deg_vals_sorted alpha_fc1_nVec];
% first sort according to the first column (the small degree value, since
% presorted above), and then sort once more according to the second degree
% value to obtain a nice ordering of the degree pairs:
degij_alpha_dbl_sorted = sortrows(degij_alpha_fc_nVec,[1 2]);
alpha_dbl_sorted       = degij_alpha_dbl_sorted(:,3);
deg_dbl_sorted         = degij_alpha_dbl_sorted(:,1:2);
boxplot(alpha_dbl_sorted,num2str(deg_dbl_sorted),'notch','on')
title(['Angle to the normal of the interface for edges with connectivity: ',num2str(1:max(deg_vals))])
xlabel('Deg of connectivity')
ylabel('Angle [deg]')

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
boxplot(alpha_dbl_sorted,num2str(deg_dbl_sorted(:,1)),'labels',M,'notch','on')
% this should be compared with the following to make sure that the counting is correct:
% h=boxplot(fc_mag_dbl_sorted,num2str(deg_dbl_sorted(:,1)));
title(['Angle to the normal of the interface for edges with connectivit: ',num2str(1:max(deg_vals))])
xlabel('Minimal deg. of connectivity')
ylabel('Angle [deg]')
clear M


% plot the stress:
figure()
% checked with part7, but the deg values should be sorted in a nicer way!
boxplot(intfStress_35,num2str(deg_vals_sorted_35),'notch','on')
title(['Interface stress for edges with connectivity: ',num2str(1:max(deg_vals))])
xlabel('Deg of connectivity')
ylabel('Interface stress [nN/um]')

figure()
% checked with part7
boxplot(lgth_vals,num2str(deg_vals_sorted),'notch','on')
title(['Interface length for edges with connectivity: ',num2str(1:max(deg_vals))])
xlabel('Deg of connectivity')
ylabel('Interface length [um]')
