% hard coded parameters:
doPrint=0;
onlyRegroup=1;



%**************************************************************************
% correlate Ecad intensity and interfacial force:
%**************************************************************************
corr_over_dl=[];
k=1;
for SIcorr_lgth=1:5:50
    goodEdgeSet=findEdges(groupedNetworks,'kPa',[8],'myoGlb',[0],'errF',500,'errs',0);
    %goodEdgeSet=findEdges(groupedNetworks,'errF',500,'errs',0);
    %[SIcorrRand_vals,Ivar_vals]=collectEdgeValues(groupedNetworks,goodEdgeSet,'SIcorr','IvarMin',2500);
    %[SIcorrRand_vals,Ivar_vals]=collectEdgeValues(groupedNetworks,goodEdgeSet,'intfRand','IvarMin',1500);
    
    %[SIcorrRand_vals]=collectEdgeValues(groupedNetworks,goodEdgeSet,'intfRand');
    %[SIcorrRand_vals]=collectEdgeValues(groupedNetworks,goodEdgeSet,'SIcorr');
    [SIcorrRand_vals,SIcorr_lgth]=collectEdgeValues(groupedNetworks,goodEdgeSet,'SIcorrRand',SIcorr_lgth,onlyRegroup);
    %[SIcorrRand_vals,SIcorr_lgth]=collectEdgeValues(groupedNetworks,goodEdgeSet,'SIcorrRand',SIcorr_lgth,onlyRegroup,'r2');

    % for getting the write dimensions:
    pixSize_mu=groupedNetworks.cluster{1}.trackedNet{7}.par.pixSize_mu;
    factor_S_to_nN_per_um=pixSize_mu*10^(-3); % See calcIntfacialStress to figure it out
    SIcorrRand_vals(:,1)= SIcorrRand_vals(:,1)*factor_S_to_nN_per_um;

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


    %calculate the numerical correlation coefficients:

    [RSIprf PSIprf RLOSIprf RUPSIprf] = corrcoef(SIcorrRand_vals(:,1),SIcorrRand_vals(:,2));
    RSIprfSTD=max(max([RSIprf-RLOSIprf RUPSIprf-RSIprf]));

    display(['corr(Fcc, Iprf): ',num2str(RSIprf(1,2),'%0.3f'),'+-',num2str(RSIprfSTD,'%0.3f'),'  (N=',num2str(length(SIcorrRand_vals(:,1))),', p=',num2str(PSIprf(1,2)),')']);
    display('!!!Since we find a significant correlation between force and intensity profile, we have achieved a subinterface force resolution!!!')
    display('The parameters of the regression are I(S):')
    display(['I = ',num2str(coeff_SI(1),'%0.1f'),' + ',num2str(coeff_SI(2)),' * S']);
    
    
    corr_over_dl(k,:)=[SIcorr_lgth RSIprf(1,2) RSIprfSTD]
    k=k+1;
end