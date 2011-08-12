function []=plotResultsForTwoStiff(groupedClusters)

goodCellSet_35=findCells(groupedClusters,'kPa',35,'myo',0,'myoGlb',[0],'errF',Inf,'errs',0);
[deg_vals_35,elE_vals_35,sumFmag_vals_35,resF_vals_35,sumFi_vals_35,sumLi_vals_35]=collectCellValues(groupedClusters,goodCellSet_35,'deg','elE','sumFmag','resF','sumFi','sumLi');

goodCellSet_8=findCells(groupedClusters,'kPa',8,'myo',0,'myoGlb',0,'errF',Inf,'errs',0);
[deg_vals_8, elE_vals_8, sumFmag_vals_8, resF_vals_8, sumFi_vals_8, sumLi_vals_8]=collectCellValues(groupedClusters,goodCellSet_8,'deg','elE','sumFmag','resF','sumFi','sumLi');

maxDeg=(max([deg_vals_35;deg_vals_8]));
G8kPa    =repmat({'8kPa'} ,length(deg_vals_8) ,1);
G35kPa   =repmat({'35kPa'},length(deg_vals_35),1);
% Gfake =repmat({''      ,maxDeg,1);

G1=vertcat(deg_vals_8,deg_vals_35);
G2=vertcat(G8kPa,G35kPa);


figure()
% checked with part7
boxplot([sumFi_vals_8;sumFi_vals_35],{G1 G2},'factorgap',[20 0],'colorgroup',G2,'color','bg','notch','on','symbol','+')
title(['Sum of interface forces for degrees: ',num2str(1:maxDeg)])
xlabel('Deg of connectivity')
ylabel('sum interface forces [nN]')
set(gca,'fontsize',16)
ylim([0,275])

figure()
% checked with part7
boxplot(sqrt(sum([resF_vals_8;resF_vals_35].^2,2)),{G1 G2},'factorgap',[20 0],'colorgroup',G2,'color','bg','notch','on','symbol','+')
title(['Residual force for degrees: ',num2str(1:maxDeg)])
xlabel('Deg of connectivity')
ylabel('Residual force [nN]')
set(gca,'fontsize',16)
ylim([0,160])

figure()
% checked with part7
boxplot([elE_vals_8;elE_vals_35],{G1 G2},'factorgap',[20 0],'colorgroup',G2,'color','bg','notch','on','symbol','+')
title(['Contraction Energy for cells with connectivity: ',num2str(1:maxDeg)])
xlabel('Deg of connectivity')
ylabel('Contraction Energy [pJ]')
set(gca,'fontsize',16)
ylim([0,0.11])


goodEdgeSet_35=findEdges(groupedClusters,'kPa',35,'myo',0,'myoGlb',[-1 0 1],'errF',Inf,'errs',0);
[deg_vals_35,lgth_vals_35,f1_vals_35,f2_vals_35,fc1_vals_35,fc2_vals_35]=collectEdgeValues(groupedClusters,goodEdgeSet_35,'deg','lgth','f1','f2','fc1','fc2');
deg_vals_sorted_35=sort(deg_vals_35,2);
fc1_mag_35 = sqrt(sum(fc1_vals_35.^2,2));

goodEdgeSet_8=findEdges(groupedClusters,'kPa',8,'myo',0,'myoGlb',0,'errF',Inf,'errs',0);
[deg_vals_8,lgth_vals_8,f1_vals_8,f2_vals_8,fc1_vals_8,fc2_vals_8]=collectEdgeValues(groupedClusters,goodEdgeSet_8,'deg','lgth','f1','f2','fc1','fc2');
deg_vals_sorted_8=sort(deg_vals_8,2);
fc1_mag_8 = sqrt(sum(fc1_vals_8.^2,2));

% figure()
degij_fcmag_35=[deg_vals_sorted_35 fc1_mag_35];
degij_fcmag_dbl_sorted_35=sortrows(degij_fcmag_35,[1 2]);
fc_mag_dbl_sorted_35 = degij_fcmag_dbl_sorted_35(:,3);
deg_dbl_sorted_35    = degij_fcmag_dbl_sorted_35(:,1:2);
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
    M{label}=['deg: ',num2str(degGroup_35),'-x; [n= ',num2str(degCount_35),']'];
    label=label+1;
end
% boxplot(fc_mag_dbl_sorted_35,num2str(deg_dbl_sorted_35(:,1)),'labels',M,'notch','on')
% title(['Interface force for edges with connectivity: ',num2str(1:max(deg_vals_35))])
% xlabel('Minimal deg. of connectivity')
% ylabel('Interface force [nN]')
% clear M


% figure()
degij_fcmag_8=[deg_vals_sorted_8 fc1_mag_8];
degij_fcmag_dbl_sorted_8=sortrows(degij_fcmag_8,[1 2]);
fc_mag_dbl_sorted_8 = degij_fcmag_dbl_sorted_8(:,3);
deg_dbl_sorted_8    = degij_fcmag_dbl_sorted_8(:,1:2);
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
    M{label}=['deg: ',num2str(degGroup_8),'-x; [n= ',num2str(degCount_8),']'];
    label=label+1;
end
% boxplot(fc_mag_dbl_sorted_8,num2str(deg_dbl_sorted_8(:,1)),'labels',M,'notch','on')
% title(['Interface force for edges with connectivity: ',num2str(1:max(deg_vals_8))])
% xlabel('Minimal deg. of connectivity')
% ylabel('Interface force [nN]')
% clear M


maxDeg=(max([deg_dbl_sorted_8(:,1);deg_dbl_sorted_35(:,1)]));
G8kPa    =repmat({'8kPa'} ,length(deg_dbl_sorted_8(:,1)) ,1);
G35kPa   =repmat({'35kPa'},length(deg_dbl_sorted_35(:,1)),1);

G1=[deg_dbl_sorted_8(:,1);deg_dbl_sorted_35(:,1)];
G2=vertcat(G8kPa,G35kPa);

figure()
% plot both together:
boxplot([fc_mag_dbl_sorted_8;fc_mag_dbl_sorted_35],{G1 G2},'factorgap',[20 0],'colorgroup',G2,'color','bg','notch','on','symbol','+')
% this should be compared with the following to make sure that the counting is correct:
% h=boxplot(fc_mag_dbl_sorted,num2str(deg_dbl_sorted(:,1)));
title(['Interface force for edges with connectivity: ',num2str(1:maxDeg)])
xlabel('Minimal deg. of connectivity')
ylabel('Interface force [nN]')
set(gca,'fontsize',16)
ylim([0,115])
clear M
