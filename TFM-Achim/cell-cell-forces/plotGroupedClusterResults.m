function []=plotGroupedClusterResults(groupedClusters)
doPrint=0;
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
% Error analysis:
%**************************************************************************
%plotErrorAnalysisNetworks(groupedClusters);

%**************************************************************************
% Estimate improvement from substracting systematic error                 :
%**************************************************************************
%plotSysError(groupedClusters)
%end %if ~onlyCorr


%**************************************************************************
% plot elastic energy and residual force over the degree.
%**************************************************************************
% goodCellSet=findCells(groupedClusters,'kPa',[8],'myo',[0],'myoGlb',[0],'errF',Inf,'errs',0);
% goodCellSet=findCells(groupedClusters,'kPa',[35],'myo',[0],'myoGlb',[-1 0 1],'type',{'ctrl','ctr','myoIIB_hp103'},'errF',Inf,'errs',0);
goodCellSet=findCells(groupedClusters,'kPa',[8],'myo',[0],'myoGlb',[0],'errF',Inf,'errs',0);
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
% checked with part7
boxplot(sumFi_vals./sumFmag_vals,deg_vals,'notch','on')
title(['Sum of interfacial forces / Sum of traction force magnitude of cells with connectivity: ',num2str(1:max(deg_vals))])
xlabel('Deg of connectivity')
ylabel('Sum of interfacial forces / Sum of traction force magnitude [1]')

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
plotResultsForTwoStiff(groupedClusters);


%**************************************************************************
% plot the interfacial force in depdence of pair degree of connectivity:
%**************************************************************************
% goodEdgeSet=findEdges(groupedClusters,'kPa',[8],'myo',[1],'type',{'myoIIA_hp93'},'errF',500,'errs',0);
% goodEdgeSet=findEdges(groupedClusters,'kPa',[8],'myo',[1],'type',{'myoIIA_hp93';'myoIIA_hp94';'myoIIB_hp103'},'errF',500,'errs',0);
goodEdgeSet=findEdges(groupedClusters,'kPa',[8],'myo',[0],'myoGlb',[0],'errF',500,'errs',0);
[deg_vals,lgth_vals,fc1_vals,nVec_vals]=collectEdgeValues(groupedClusters,goodEdgeSet,'deg','lgth','fc1','nVec');
deg_vals_sorted=sort(deg_vals,2);
fc1_mag    = sqrt(sum(fc1_vals.^2,2));
intfStress = fc1_mag./lgth_vals;

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
boxplot(intfStress,num2str(deg_vals_sorted),'notch','on')
title(['Interface stress for edges with connectivity: ',num2str(1:max(deg_vals))])
xlabel('Deg of connectivity')
ylabel('Interface stress [nN/um]')

figure()
% checked with part7
boxplot(lgth_vals,num2str(deg_vals_sorted),'notch','on')
title(['Interface length for edges with connectivity: ',num2str(1:max(deg_vals))])
xlabel('Deg of connectivity')
ylabel('Interface length [um]')


figure()
% checked with part7
degij_stress=[deg_vals_sorted intfStress];
% first sort according to the first column (the small degree value, since
% presorted above), and then sort once more according to the second degree
% value to obtain a nice ordering of the degree pairs:
degij_stress_dbl_sorted=sortrows(degij_stress,[1 2]);
stress_dbl_sorted = degij_stress_dbl_sorted(:,3);
deg_dbl_sorted    = degij_stress_dbl_sorted(:,1:2);
boxplot(stress_dbl_sorted,num2str(deg_dbl_sorted),'notch','on')
title(['Interface stress for edges with connectivity: ',num2str(1:max(deg_vals))])
xlabel('Deg of connectivity')
ylabel('Interface stress [nN]')

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
boxplot(stress_dbl_sorted,num2str(deg_dbl_sorted(:,1)),'labels',M,'notch','on')
% this should be compared with the following to make sure that the counting is correct:
% h=boxplot(fc_mag_dbl_sorted,num2str(deg_dbl_sorted(:,1)));
title(['Interface stress for edges with connectivity: ',num2str(1:max(deg_vals))])
xlabel('Minimal deg. of connectivity')
ylabel('Interface stress [nN/um]')
clear M

%**************************************************************************
% plot the angle in depdence of pair degree of connectivity:
%**************************************************************************
% goodEdgeSet=findEdges(groupedClusters,'kPa',[8],'myo',[1],'type',{'myoIIA_hp93'},'errF',500,'errs',0);
% goodEdgeSet=findEdges(groupedClusters,'kPa',[8],'myo',[1],'type',{'myoIIA_hp93';'myoIIA_hp94';'myoIIB_hp103'},'errF',500,'errs',0);
%goodEdgeSet=findEdges(groupedClusters,'kPa',[8],'myo',[0],'myoGlb',[0],'errF',500,'errs',0);
% goodEdgeSet=findEdges(groupedClusters,'kPa',[8],'errF',500,'errs',0);
goodEdgeSet=findEdges(groupedClusters,'kPa',[8],'myo',[0],'myoGlb',[0],'errF',30,'errs',0);
[deg_vals,lgth_vals,fc1_vals,nVec_vals]=collectEdgeValues(groupedClusters,goodEdgeSet,'deg','lgth','fc1','nVec');
deg_vals_sorted=sort(deg_vals,2);
fc1_mag = sqrt(sum(fc1_vals.^2,2));
alpha_fc1_nVec= acosd(dot(fc1_vals,nVec_vals,2)./(fc1_mag));

%quick cheat:
cv=alpha_fc1_nVec>90;
alpha_fc1_nVec(cv)=180-alpha_fc1_nVec(cv);

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

%end %if ~onlyCorr
%**************************************************************************
% Plot ctr-ctr, ctr-myo, myo-myo, for degree 1-1 interfaces for E=8,35kPa:*
%**************************************************************************
plotIntForceDeg11(groupedClusters);

%end %if ~onlyCorr

%**************************************************************************
% correlate Ecad intensity and interfacial force:
%**************************************************************************
% ech single conditions works really well, but mixtures are a bit worse.
% goodEdgeSet=findEdges(groupedClusters,'kPa',[8],'myoGlb',[0],'errF',500,'errs',0);
% goodEdgeSet=findEdges(groupedClusters,'kPa',[8],'myoGlb',[1],'errF',500,'errs',0);
goodEdgeSet=findEdges(groupedClusters,'kPa',[8],'myoGlb',[0],'errF',500,'errs',0);
% goodEdgeSet=findEdges(groupedClusters,'kPa',[8],'errF',500,'errs',0);

[fc1_vals,Itot_vals,Iavg_vals,SIcorr_vals]=collectEdgeValues(groupedClusters,goodEdgeSet,'fc1','Itot','Iavg','SIcorr');
fc1_mag = sqrt(sum(fc1_vals.^2,2));
% bin the data
edgesFBins=linspace(0,max(fc1_mag),20);
[n,bin]=histc(fc1_mag,edgesFBins);
corrFItot=[];
for iBin=1:max(bin)
    checkVec= (iBin==bin);
    corrFItot(iBin).FVals    =fc1_mag(checkVec);
    corrFItot(iBin).FValsMean=mean(corrFItot(iBin).FVals);
    corrFItot(iBin).FValsSTD = std(corrFItot(iBin).FVals);
    
    corrFItot(iBin).IVals    = Itot_vals(checkVec);
    corrFItot(iBin).IValsMean=mean(corrFItot(iBin).IVals);
    corrFItot(iBin).IValsSTD = std(corrFItot(iBin).IVals);
end

figure()
plot(fc1_mag,Itot_vals,'.k','MarkerSize',6)
hold on;
errorbarxy([corrFItot.FValsMean],[corrFItot.IValsMean],[corrFItot.FValsSTD],[corrFItot.IValsSTD],[],[],'sk','k');
hold on;
% plot robust fit:
[fc1_mag id] =sortrows(fc1_mag);
Itot_vals = Itot_vals(id);
coeff_FI = regress(Itot_vals,[ones(size(Itot_vals)) fc1_mag]);
%coeff = robustfit(fc1_mag,Itot_vals);
plot(fc1_mag,coeff_FI(1)+coeff_FI(2)*fc1_mag,'k','LineWidth',2)
title('Correlation Ecad intensity / cell-cell forces')
xlabel('Cell-cell force magnitude [nN]')
ylabel('Integrated Ecad intensity [a.u.]')
box on
set(gca,'LineWidth',2,'FontSize',20)
hold off;
if doPrint
    saveas(gcf,'corrFI.fig');
    saveas(gcf,'corrFI.eps','psc2');
end


% figure()
% plot(fc1_mag,Iavg_vals,'o')
% title('Correlation Ecad intensity / cell-cell forces')
% xlabel('Cell-cell force magnitude [nN]')
% ylabel('Average Ecad intensity [a.u.]')

pixSize_mu=groupedClusters.cluster{1}.trackedNet{1}.par.pixSize_mu;
factor_S_to_nN_per_um=pixSize_mu*10^(-3); % See calcIntfacialStress to figure it out
SIcorr_vals(:,1)=SIcorr_vals(:,1)*factor_S_to_nN_per_um;

edgesSBins=linspace(0,max(SIcorr_vals(:,1)),20);
[n,bin]=histc(SIcorr_vals(:,1),edgesSBins);
corrSItot=[];
for iBin=1:max(bin)
    checkVec= (iBin==bin);
    corrSItot(iBin).SVals    =SIcorr_vals(checkVec,1);
    corrSItot(iBin).SValsMean=mean(corrSItot(iBin).SVals);
    corrSItot(iBin).SValsSTD = std(corrSItot(iBin).SVals);
    
    corrSItot(iBin).IVals    =SIcorr_vals(checkVec,2);
    corrSItot(iBin).IValsMean=mean(corrSItot(iBin).IVals);
    corrSItot(iBin).IValsSTD = std(corrSItot(iBin).IVals);
end

figure()
plot(SIcorr_vals(:,1),SIcorr_vals(:,2),'.k','MarkerSize',3)
hold on;
errorbarxy([corrSItot.SValsMean],[corrSItot.IValsMean],[corrSItot.SValsSTD],[corrSItot.IValsSTD],[],[],'sk','k');
hold on;
% plot robust fit:
SIcorr_vals =sortrows(SIcorr_vals);
coeff_SI = regress(SIcorr_vals(:,2),[ones(size(SIcorr_vals(:,2))) SIcorr_vals(:,1)]);
% coeff = robustfit(SIcorr_vals(:,1),SIcorr_vals(:,2));
plot(SIcorr_vals(:,1),coeff_SI(1)+coeff_SI(2)*SIcorr_vals(:,1),'k','LineWidth',2)
title('Correlation Ecad intensity profile/ cell-cell stress profile')
xlabel('Cell-cell stress [nN/um]')
ylabel('Locally average Ecad intensity [a.u.]')
box on
set(gca,'LineWidth',2,'FontSize',20)
hold off;
if doPrint
    saveas(gcf,'corrSI.fig');
    saveas(gcf,'corrSI.eps','psc2');
end


%calculate the numerical correlation coefficients:
[RFItot PFItot RLOItot RUPItot] = corrcoef(fc1_mag,Itot_vals);
RFItotSTD=max(max([RFItot-RLOItot RUPItot-RFItot]));

[RFIavg PFIavg RLOIavg RUPIavg] = corrcoef(fc1_mag,Iavg_vals);
RFIavgSTD=max(max([RFIavg-RLOIavg RUPIavg-RFIavg]));

[RSIprf PSIprf RLOSIprf RUPSIprf] = corrcoef(SIcorr_vals(:,1),SIcorr_vals(:,2));
RSIprfSTD=max(max([RSIprf-RLOSIprf RUPSIprf-RSIprf]));

display(['corr(Fcc, Itot): ',num2str(RFItot(1,2),'%0.3f'),'+-',num2str(RFItotSTD,'%0.3f'),'  (N=',num2str(length(fc1_mag),'% 6.0f')         ,', p=',num2str(PFItot(1,2)),')']);
display(['corr(Fcc, Iavg): ',num2str(RFIavg(1,2),'%0.3f'),'+-',num2str(RFIavgSTD,'%0.3f'),'  (N=',num2str(length(fc1_mag))         ,', p=',num2str(PFIavg(1,2)),')']);
display(['corr(Fcc, Iprf): ',num2str(RSIprf(1,2),'%0.3f'),'+-',num2str(RSIprfSTD,'%0.3f'),'  (N=',num2str(length(SIcorr_vals(:,1))),', p=',num2str(PSIprf(1,2)),')']);
display('!!!Since we find a significant correlation between force and intensity profile, we have achieved a subinterface force resolution!!!')
display('The parameters of the regression are I(F):')
display(['I = ',num2str(coeff_FI(1),'%0.1f'),' + ',num2str(coeff_FI(2),'%0.1f'),' * F']);
display('The parameters of the regression are I(S):')
display(['I = ',num2str(coeff_SI(1),'%0.1f'),' + ',num2str(coeff_SI(2)),' * S']);

%end %if ~onlyCorr

normVar=1;
tBtwFrms=240;
aveType='nanmean'; % first checks 'none', 'nanmean', 'mean' makes little difference
maxLag =round(7200/tBtwFrms); % round(7200/tBtwFrms) means a maxLag of 2h
relErrF_val_corr=Inf;
%**************************************************************************
% correlate forces for control cells:
%**************************************************************************
% goodCellSet   = findCells(groupedClusters,'kPa',35,'deg',[2 3 4 5 6 7],'myo',[0],'divGlb',[-1 0 1],'relErrF',relErrF_val_corr,'errs',0);
goodCellSet   = findCells(groupedClusters,'kPa',[8],'deg',[2 3 4 5 6 7],'myo',[0],'myoGlb',[-1 0 1],'relErrF',relErrF_val_corr,'errs',0);
[corrSets]    = collectCellValues(groupedClusters,goodCellSet,'corr');
[corrResults] = calCorrResults(corrSets,maxLag,'usefm',normVar,tBtwFrms,aveType);


%**************************************************************************
% correlate forces for myosin cells:
%**************************************************************************
% goodCellSet=findCells(groupedClusters,'kPa',[8],'deg',[2 3 4 5 6 7],'myo',[1],'type',{'tln1'},'relErrF',relErrF_val_corr,'errs',0);
% goodCellSet=findCells(groupedClusters,'kPa',35,'deg',[2 3 4 5 6 7],'myo',[1],'type',{'myoIIB_hp103'},'errF',errF_val_corr,'errs',0);
goodCellSet=findCells(groupedClusters,'kPa',[8],'deg',[2 3 4 5 6 7],'myo',[1],'divGlb',[-1 0 1],'type',{'myoIIA_hp93';'myoIIA_hp94'},'relErrF',relErrF_val_corr,'errs',0);
if ~isempty(goodCellSet) && ~isempty(goodCellSet(1).cellId)
    [corrSets]=collectCellValues(groupedClusters,goodCellSet,'corr');
    [corrResults]=calCorrResults(corrSets,maxLag,'usefm',normVar,tBtwFrms,aveType);
else
    display('No myosin cells of this type found!')
end

end %if ~onlyCorr
%**************************************************************************
% correlate forces and Ecad intensity:
%**************************************************************************

normVar=1;
tBtwFrms=150;
aveType='nanmean'; % first checks 'none', 'nanmean', 'mean' makes little difference
% plots generated for Gordon Conference used 'none' but 'nanmean' gives
% ~equal results
maxLag =round(3600/tBtwFrms); % round(7200/tBtwFrms) means a maxLag of 2h
relErrF_val_corr=Inf;

%goodCellSet=findCells(groupedClusters,'kPa',[8],'deg',[2 3 4 5 6 7],'myo',[1],'type',{'tln1'},'errF',errF_val_corr,'errs',0);
%goodCellSet=findCells(groupedClusters,'kPa',35,'deg',[2 3 4 5 6 7],'myo',[1],'type',{'myoIIB_hp103'},'errF',errF_val_corr,'errs',0);
%goodEdgeSet=findEdges(groupedClusters,'kPa',[8],'asmbly',[1],'relErrF',relErrF_val_corr,'errs',0);
% goodEdgeSet=findEdges(groupedClusters,'kPa',[8],'asmbly',[1],'relErrF',relErrF_val_corr,'errs',0);
% goodEdgeSet=findEdges(groupedClusters,'kPa',[8],'asmbly',[1],'relErrF',relErrF_val_corr,'errs',0);
% goodEdgeSet=findEdges(groupedClusters,'myo',[1],'relErrF',relErrF_val_corr,'errs',0);
goodEdgeSet=findEdges(groupedClusters,'myo',[0],'type',{'myoIIA_hp93';'myoIIA_hp94'},'relErrF',relErrF_val_corr,'errs',0);
%goodEdgeSet=findEdges(groupedClusters,'relErrF',relErrF_val_corr,'errs',0,'dt',280);
if ~isempty(goodEdgeSet) && ~isempty(goodEdgeSet(1).edgeId)
    [corrSets]=collectEdgeValues(groupedClusters,goodEdgeSet,'corr');
    %[corrResults]=calCorrResultsInt(corrSets,maxLag,'usefm',normVar,tBtwFrms,aveType,'useItot');
    [corrResults]=calCorrResultsInt(corrSets,maxLag,'usefm',normVar,tBtwFrms,aveType,'useItot');
else
    display('No myosin cells of this type found!')
end

%**************************************************************************
% Plot the edges from above:
%**************************************************************************
%goodEdgeSet=findEdges(groupedClusters,'kPa',[8],'asmbly',[-1],'relErrF',relErrF_val_corr,'errs',0);