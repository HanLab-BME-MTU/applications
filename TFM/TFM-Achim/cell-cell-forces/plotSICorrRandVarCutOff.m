% hard coded parameters:
doPrint=0;
onlyRegroup=0;



%**************************************************************************
% correlate Ecad intensity and interfacial force:
%**************************************************************************
corr_over_IvarMin=[];
k=1;
for IvarMinVal=0:100:5000
    goodEdgeSet=findEdges(groupedNetworks,'kPa',[8],'myoGlb',[0],'errF',500,'errs',0);
    [SIcorr_vals    ,Ivar_vals]=collectEdgeValues(groupedNetworks,goodEdgeSet,'SIcorr'  ,'IvarMin',IvarMinVal);
    [SIcorrRand_vals,Ivar_vals]=collectEdgeValues(groupedNetworks,goodEdgeSet,'intfRand','IvarMin',IvarMinVal);

    % for getting the write dimensions:
    pixSize_mu=groupedNetworks.cluster{1}.trackedNet{7}.par.pixSize_mu;
    factor_S_to_nN_per_um=pixSize_mu*10^(-3); % See calcIntfacialStress to figure it out
    SIcorr_vals(:,1)    = SIcorr_vals(:,1)*factor_S_to_nN_per_um;
    SIcorrRand_vals(:,1)= SIcorrRand_vals(:,1)*factor_S_to_nN_per_um;

    % This is only needed to plot the errorbarxy
    edgesSBins=linspace(0,max(SIcorrRand_vals(:,1)),20);
    [n,bin]=histc(SIcorrRand_vals(:,1),edgesSBins);
    corrSItot=[];
    for iBin=1:max(bin)
        checkVec= (iBin==bin);
        corrSItot(iBin).SVals    =SIcorrRand_vals(checkVec,1);
        corrSItot(iBin).SValsMean=mean(corrSItot(iBin).SVals);
        corrSItot(iBin).SValsSTD = std(corrSItot(iBin).SVals);

        corrSItot(iBin).IVals    =SIcorrRand_vals(checkVec,2);
        corrSItot(iBin).IValsMean=mean(corrSItot(iBin).IVals);
        corrSItot(iBin).IValsSTD = std(corrSItot(iBin).IVals);
    end

    figure()
    plot(SIcorrRand_vals(:,1),SIcorrRand_vals(:,2),'.k','MarkerSize',3)
    hold on;
    errorbarxy([corrSItot.SValsMean],[corrSItot.IValsMean],[corrSItot.SValsSTD],[corrSItot.IValsSTD],[],[],'sk','k');
    hold on;
    % plot robust fit:
    SIcorrRand_vals =sortrows(SIcorrRand_vals);
    coeff_SI = regress(SIcorrRand_vals(:,2),[ones(size(SIcorrRand_vals(:,2))) SIcorrRand_vals(:,1)]);
    % coeff = robustfit(SIcorrRand_vals(:,1),SIcorrRand_vals(:,2));
    plot(SIcorrRand_vals(:,1),coeff_SI(1)+coeff_SI(2)*SIcorrRand_vals(:,1),'k','LineWidth',2)
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

    %calculate the numerical correlation coefficients control:
    [RSIprf_ctr PSIprf_ctr RLOSIprf_ctr RUPSIprf_ctr] = corrcoef(SIcorr_vals(:,1),SIcorr_vals(:,2));
    RSIprf_ctrSTD=max(max([RSIprf_ctr-RLOSIprf_ctr RUPSIprf_ctr-RSIprf_ctr]));

    display(['corr(Fcc, Iprf): ',num2str(RSIprf_ctr(1,2),'%0.3f'),'+-',num2str(RSIprf_ctrSTD,'%0.3f'),'  (N=',num2str(length(SIcorr_vals(:,1))),', p=',num2str(PSIprf_ctr(1,2)),')']);
    display('!!!Since we find a significant correlation between force and intensity profile, we have achieved a subinterface force resolution!!!')

    %calculate the numerical correlation coefficients randomized:
    [RSIprf PSIprf RLOSIprf RUPSIprf] = corrcoef(SIcorrRand_vals(:,1),SIcorrRand_vals(:,2));
    RSIprfSTD=max(max([RSIprf-RLOSIprf RUPSIprf-RSIprf]));

    display(['corr(Fcc, Iprf): ',num2str(RSIprf(1,2),'%0.3f'),'+-',num2str(RSIprfSTD,'%0.3f'),'  (N=',num2str(length(SIcorrRand_vals(:,1))),', p=',num2str(PSIprf(1,2)),')']);
    display('!!!Since we find a significant correlation between force and intensity profile, we have achieved a subinterface force resolution!!!')
    display('The parameters of the regression are I(S):')
    display(['I = ',num2str(coeff_SI(1),'%0.1f'),' + ',num2str(coeff_SI(2)),' * S']);
    
    
    corr_over_IvarMin(k,:)=[IvarMinVal RSIprf_ctr(1,2) RSIprf_ctrSTD  RSIprf(1,2) RSIprfSTD]
    k=k+1;
end

figure()
errorbar(corr_over_IvarMin(:,1),corr_over_IvarMin(:,2),corr_over_IvarMin(:,3),'-k')
hold on
errorbar(corr_over_IvarMin(:,1),corr_over_IvarMin(:,4),corr_over_IvarMin(:,5),'-r')
plot(corr_over_IvarMin(:,1),corr_over_IvarMin(:,4)./corr_over_IvarMin(:,2),'-b')