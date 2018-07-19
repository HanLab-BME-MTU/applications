% hard coded parameters:
kPa=8;
dpts=5; % usual sampling dpts=5; for dense sampling use dpts=2
%SIcorr_lgth_plot_list=[1, 5, 10, 20, 35];
SIcorr_lgth_plot_list=[1, 6, 11, 21, 36];

bin_vec=linspace(-1,1,21); % density of histogram bins (resolution of correlation coefficient)
marker=[{'-k'},{'-ro'},{'-bs'},{'-m*'},{'-c+'},{'-gd'},{'-yo'},{'-ks'},{'-.ro'},{'-.bs'},{'-.m*'},{'-.c+'},{'-.gd'},{'-.yo'},{'-.ks'}];
marker=horzcat(marker,[{'--k'},{'--ro'},{'--bs'},{'--m*'},{'--c+'},{'--gd'},{'--yo'},{'--ks'},{':ro'},{':bs'},{':m*'},{':c+'},{':gd'},{':yo'},{':ks'}]);

% do not touch:
rng('default'); % fixed seed so the randomization give always the same value!

corr_of_index_rand   =[];
corr_of_index_correct=[];
h=[];
M=[];
k=1;
close all;

% only extract the scale for dl in um:
%goodEdgeSet=findEdges(groupedNetworks,'kPa',kPa,'myoGlb',[0],'errF',Inf,'errs',0); % ('errF',500) gives exactly the same results
goodEdgeSet=findEdges(groupedNetworks,'myoGlb',[0],'errF',Inf,'errs',0);
[dlpix_vals]=collectEdgeValues(groupedNetworks,goodEdgeSet,'dlpix');
% hist(dlpix_vals,100)
dl_um=mean(dlpix_vals)*groupedNetworks.cluster{1}.trackedNet{1}.par.pixSize_mu;
% median(dlpix_vals)

%**************************************************************************
% correlate Ecad intensity and interfacial force after full randomization
%**************************************************************************
% plot the full randomization first:
goodEdgeSet=findEdges(groupedNetworks,'kPa',kPa,'myoGlb',[0],'errF',Inf,'errs',0); % ('errF',500) gives exactly the same results
[SIcorrCoef_vals_intfRand]=collectEdgeValues(groupedNetworks,goodEdgeSet,'SIcorrCoefIntfRand');

% prepare the normalized histogram
num_pairs_intfRand=length(SIcorrCoef_vals_intfRand(:,1));
n_in_bins_intfRand=hist(SIcorrCoef_vals_intfRand(:,1),bin_vec);

    
for onlyRegroup=0:1
    % for plotting:
    h=[];
    M=[];    
    k=1;
    markerId=1;
    
    % do some statistics first (the last entries #6 and #7 hold the
    % clusterIds and edgeIds, respectively):
    % 1) count the number of clusters that were used
    clusterList_intfRand=unique(SIcorrCoef_vals_intfRand(:,6));
    num_clusters_intfRand=length(clusterList_intfRand);
    % 1) count the number of edges that were used
    edgeList_intfRand=unique(SIcorrCoef_vals_intfRand(:,6:7),'rows');
    num_edges_intfRand=length(edgeList_intfRand(:,2));
    
    figure(onlyRegroup+1)
    % plot the normalized histogram
    currh=plot(bin_vec,n_in_bins_intfRand/num_pairs_intfRand,marker{markerId},'LineWidth',4);
    hold on;
    h=horzcat(h,currh(1)); % 1 for marker, 2 for colored line;
    M{markerId}=['dl= full;  ',';  median= ',num2str(median(SIcorrCoef_vals_intfRand(:,1)),'%.3f'),';  #vals= ', num2str(num_pairs_intfRand),';  #intfs= ', num2str(num_edges_intfRand),';  #clusts= ', num2str(num_clusters_intfRand)];

    if onlyRegroup
        corr_of_index_correct(k,:)=[NaN num_pairs_intfRand median(SIcorrCoef_vals_intfRand(:,1)) mean(SIcorrCoef_vals_intfRand(:,1)) num_clusters_intfRand num_edges_intfRand];
    else
        corr_of_index_rand(k,:)=[NaN num_pairs_intfRand median(SIcorrCoef_vals_intfRand(:,1)) mean(SIcorrCoef_vals_intfRand(:,1)) num_clusters_intfRand num_edges_intfRand];
    end
    k=k+1;
    markerId=markerId+1;
    
    %**************************************************************************
    % correlate Ecad intensity and interfacial force:
    %**************************************************************************
    for SIcorr_lgth=1:dpts:50   %50 is enough, there are no longer interfaces
        goodEdgeSet=findEdges(groupedNetworks,'kPa',kPa,'myoGlb',[0],'errF',Inf,'errs',0); % ('errF',500) gives exactly the same results    
        [SIcorrCoef_vals,~,~]=collectEdgeValues(groupedNetworks,goodEdgeSet,'SIcorrCoefRand',SIcorr_lgth,onlyRegroup);

        if ~isempty(SIcorrCoef_vals)
            % do some statistics first (the last entries #6 and #7 hold the
            % clusterIds and edgeIds, respectively):
            % 1) count the number of clusters that were used
            clusterList=unique(SIcorrCoef_vals(:,6));
            num_clusters=length(clusterList);
            % 1) count the number of edges that were used
            edgeList=unique(SIcorrCoef_vals(:,6:7),'rows');
            num_edges=length(edgeList(:,2));
            
            % do only plot the histograms we want to plot (not too many otherwise it's getting too confusing)
            if find(SIcorr_lgth==SIcorr_lgth_plot_list)
                % prepare the normalized histogram
                num_pairs=length(SIcorrCoef_vals(:,1));
                n_in_bins=hist(SIcorrCoef_vals(:,1),bin_vec);

                % plot the normalized histogram
                currh=plot(bin_vec,n_in_bins/num_pairs,marker{markerId},'LineWidth',2);
                h=horzcat(h,currh(1)); % 1 for marker, 2 for colored line;
                hold on;
                xlabel('correlation coefficient')
                ylabel('probability')
                box on
                set(gca,'LineWidth',2,'FontSize',20)
                M{markerId}=['dl= ',num2str(SIcorr_lgth),';  median= ',num2str(median(SIcorrCoef_vals(:,1)),'%.3f'),';  #vals= ', num2str(num_pairs),';  #intfs= ', num2str(num_edges),';  #clusts= ', num2str(num_clusters)];
                legend(h,M);
                if onlyRegroup
                    title('Distribution of corr-coefs (correct pairings)')
                else
                    title('Distribution of corr-coefs (randomized)')
                end
                
                markerId=markerId+1;
            end  
            
            if onlyRegroup
                corr_of_index_correct(k,:)=[SIcorr_lgth num_pairs median(SIcorrCoef_vals(:,1)) mean(SIcorrCoef_vals(:,1)) num_clusters num_edges];
            else
                corr_of_index_rand(k,:)=[SIcorr_lgth num_pairs median(SIcorrCoef_vals(:,1)) mean(SIcorrCoef_vals(:,1)) num_clusters num_edges];
            end
            k=k+1;
        end
    end
end
   
figure(3)
h=[];
M=[]; 
max_dl=dl_um*max(corr_of_index_rand(:,1));
plot([-max_dl,max_dl],[0.5,0.5],'-k')
hold on
% the first entry is the full randomization so we start at the second 
% entry (2:end):
currh=plot(dl_um*corr_of_index_correct(2:end,1),corr_of_index_correct(2:end,3),'-ok','LineWidth',2);
h=horzcat(h,currh(1));
M{1}=['correct pairings'];
currh=plot(dl_um*corr_of_index_rand(2:end,1),corr_of_index_rand(2:end,3),'-ob','LineWidth',2);
h=horzcat(h,currh(1));
M{2}=['randomized pairings'];
currh=plot(dl_um*corr_of_index_rand(2:end,1),corr_of_index_rand(2:end,3)./corr_of_index_correct(2:end,3),'-or','LineWidth',2);
h=horzcat(h,currh(1));
M{3}=['ratio: random/correct'];
set(gca,'LineWidth',2,'FontSize',20)
xlabel(['length of randomization interval [um]'])
%xlabel(['length of randomization interval [',num2str(dl_um,'%.2f'),'um]'])
ylabel('median correlation coefficient')
box on
legend(h,M);
axis([0,dl_um*corr_of_index_correct(end,1),-0.2,1.1])


% print the histogram of the interface length:
figure(4)
goodEdgeSet=findEdges(groupedNetworks,'kPa',kPa,'myoGlb',[0],'errF',Inf,'errs',0); % ('errF',500) gives exactly the same results
[legth_vals]=collectEdgeValues(groupedNetworks,goodEdgeSet,'lgth');
hist(legth_vals,20)
set(gca,'LineWidth',2,'FontSize',20)
xlabel(['interface length [um]'])
ylabel(['counts'])
box on
xmax=63;
axis([0,xmax,0.0,300])
text(xmax-25,275,['median= ',num2str(median(legth_vals),'%.1f')],'FontSize',20)
text(xmax-25,250,['mean= ',num2str(mean(legth_vals),'%.1f')],'FontSize',20)
title(['distribution of interface length: ',num2str(kPa),'kPa'])


%display the results as table:
display('correct pairings:')
num2str(corr_of_index_correct(2:end,:))
display('randomized pairings:')
num2str(corr_of_index_rand(2:end,:))
display('full randomization:')
num2str(corr_of_index_rand(1,:))