
%%
All_ChMP_feature_thisML_ch1 = nan(33,8,384);
All_ChMP_feature_thisML_ch2 = nan(33,8,384);

% get the control medians    
for iMD = 1 : 384
    All_ChMP_feature_thisML_ch1(:,:,iMD) = ChMP_feature_thisML{iMD}{1}(:,:);
    All_ChMP_feature_thisML_ch2(:,:,iMD) = ChMP_feature_thisML{iMD}{2}(:,:);
end

 all_ch1_median = squeeze(All_ChMP_feature_thisML_ch1(:,5,:));
 all_ch1_mean = squeeze(All_ChMP_feature_thisML_ch1(:,1,:));
 all_ch2_median = squeeze(All_ChMP_feature_thisML_ch2(:,5,:));
 all_ch2_mean = squeeze(All_ChMP_feature_thisML_ch2(:,1,:));
 
 for iF = 1 : 8
    all_ch1_status{iF} = squeeze(All_ChMP_feature_thisML_ch1(:,iF,:));
    all_ch2_status{iF} = squeeze(All_ChMP_feature_thisML_ch2(:,iF,:));
    
    all_ch1_median_valid = all_ch1_median(find(feature_index>0),:);
 
 end
 
 
 all_ch1_median_valid = all_ch1_median(find(feature_index>0),:);
 all_ch1_mean_valid = all_ch1_mean(find(feature_index>0),:);
 all_ch2_median_valid =all_ch2_median(find(feature_index>0),:);
 all_ch2_mean_valid =all_ch2_mean(find(feature_index>0),:);
 
 
 for i = 1 : 12
     
    all_ch1_median_valid_normalized(i,:) = (all_ch1_median_valid(i,:)...
        -nanmean(all_ch1_median_valid(i,:)))...
        ./nanstd(all_ch1_median_valid(i,:));
     all_ch1_mean_valid_normalized(i,:) = (all_ch1_mean_valid(i,:)...
        -nanmean(all_ch1_mean_valid(i,:)))...
        ./nanstd(all_ch1_median_valid(i,:));
      all_ch2_median_valid_normalized(i,:) = (all_ch2_median_valid(i,:)...
        -nanmean(all_ch2_median_valid(i,:)))...
        ./nanstd(all_ch1_median_valid(i,:));
      all_ch2_mean_valid_normalized(i,:) = (all_ch2_mean_valid(i,:)...
        -nanmean(all_ch2_mean_valid(i,:)))...
        ./nanstd(all_ch2_mean_valid(i,:));
   
     
 end

 all_meadian_matrix = [all_ch1_median_valid_normalized; all_ch2_median_valid_normalized];
 
 save('all_ch_1_2_true_median_mean.mat');
 
 
 A_IDs = [24*(1:14)+1 3*24+10];
 N_IDs = [24*(1:14)+24 9*24+15];
 C_IDs = 1:384;
 C_IDs = setdiff(setdiff(C_IDs,A_IDs),N_IDs);
 
 
 figure(1);
 hold off; plot3(all_ch1_median_valid_normalized(2,A_IDs),all_ch1_median_valid_normalized(9,A_IDs),all_ch2_median_valid_normalized(2,A_IDs),'r.');
 hold on; plot3(all_ch1_median_valid_normalized(2,N_IDs),all_ch1_median_valid_normalized(9,N_IDs),all_ch2_median_valid_normalized(2,N_IDs),'g.');
 hold on; plot3(all_ch1_median_valid_normalized(2,C_IDs),all_ch1_median_valid_normalized(9,C_IDs),all_ch2_median_valid_normalized(2,C_IDs),'b.');
 xlabel('vim length');
 ylabel('vim st');
 zlabel('mt length');
 
 
 [pc,score,latent,tsquare]=princomp(all_meadian_matrix');
 
 pced_data = (all_meadian_matrix'*pc)';
 
 figure(2);
 hold off; plot3(pced_data(1,A_IDs),pced_data(2,A_IDs),pced_data(3,A_IDs),'r.');
 hold on; plot3(pced_data(1,N_IDs),pced_data(2,N_IDs),pced_data(3,N_IDs),'g.');
 hold on; plot3(pced_data(1,C_IDs),pced_data(2,C_IDs),pced_data(3,C_IDs),'b.');
 xlabel('pc1');
 ylabel('pc2');
 zlabel('pc3');
 
 a=1;
 
 
 
 
 
%    
% 
% for iF = 1:32
%     if iF==2    
%     else
%         
%     control_median = Control_ChMP_feature_thisML_ch1(iF,5,valid_control_ID);
%     end
%     CFMP_feature_50p_wells(iF) = prctile(control_median(:),50);
%     
%     CFMP_feature_5p_wells(iF) = prctile(control_median(:),5);
%     
%     CFMP_feature_95p_wells(iF) = prctile(control_median(:),95);
%     
%     CFMP_feature_2p_wells(iF) = prctile(control_median(:),2);
%     
%     CFMP_feature_98p_wells(iF) = prctile(control_median(:),98);    
%     
%     CFMP_feature_25p_wells(iF) = prctile(control_median(:),25);
%     
%     CFMP_feature_75p_wells(iF) = prctile(control_median(:),75);
% end
%   
% 
% textcell = {'Straightness of filament (for each filament)',... % 1
%     'Length (for each filament)',...                           % 2
%     'Pixel number of segmentation (for each filament)',...     % 3
%     'Filament density (for each valid pixel)',...                   % 4
%     'Scrabled filament density (for each valid pixel)',...          % 5
%     'Filament orientation (for each pixel)',...                % 6
%     'Filament orientation (for each pixel,centered)',...       % 7
%     'Filamet intensity (integrated for each filament)',...           % 8
%     'Filamet intensity (average for each filament)',...              % 9
%     'Filamet intensity (integrated for each dilated filament)',...   % 10
%     'Filamet intensity (average for each dilated filament)',...      % 11
%     'Scale of detected filaments (for each pixel)',...                   % 12
%     'ST Responce (integrated for each filament)',...            % 13
%     'ST Responce (average for each filament)',...               % 14
%     'ST Responce (integrated for each dilated filament)',...    % 15
%     'ST Responce (average for each dilated filament)',...       % 16
%     'Filament curvature (average for each filament)',...     % 17
%     'Filament curvature (for each pixel on filaments)',...   % 18
%     'NMS Sum Periphery-Center Ratio (for each sections for cell (Dpc40p))',...           % 19
%     'NMS Ave Periphery-Center Ratio (for each sections for cell(Dpc40p))',...            % 20
%     'Intensity Sum Periphery-Center Ratio (for each sections for cell(Dpc40p))',...      % 21
%     'Intensity Ave Periphery-Center Ratio (for each sections for cell(Dpc40p))',...      % 22
%     'Filament Pixel Sum Periphery-Center Ratio (for each sections for cell(Dpc40p))',... % 23
%     'Filament Density Periphery-Center Ratio (for each sections for cell(Dpc40p))',...   % 24
%     'NMS Sum Periphery-Center Ratio (for each sections for cell (AutoDpc))',...           % 25
%     'NMS Ave Periphery-Center Ratio (for each sections for cell(AutoDpc))',...            % 26
%     'Intensity Sum Periphery-Center Ratio (for each sections for cell(AutoDpc))',...      % 27
%     'Intensity Ave Periphery-Center Ratio (for each sections for cell(AutoDpc))',...      % 28
%     'Filament Pixel Sum Periphery-Center Ratio (for each sections for cell(AutoDpc))',... % 29
%     'Filament Density Periphery-Center Ratio (for each sections for cell(AutoDpc))',...   % 30
%     'Filament Centripetal angle (for each filament)',...   % 31
%     'Filament Centripetal angle (for each pixel)',...      % 32
%     'Number of Nucleus(per well)'...                    %33
%     };
% 
% Group_ROOT_DIR='/project/cellbiology/gdanuser/vimentin/ding/fromTony/Control_P567/MovieList_P567';

% nMovie = 200+30+10;
% Group_ROOT_DIR='.';
all_hits_pool = [];
hits_count = 0;
control_hits_count = 0;
        control_hits_pool=[];
        
for iF =  [ 1:3 6:17 19:32] 
    % for the background, plot 25%~50%~75% patch
    % then the mean
    h1 = figure(1); close;
    h1 = figure(1); hold off;
    
    p1 = patch([1  numel(valid_control_ID) numel(valid_control_ID) 1 1],...
        [CFMP_feature_2p_wells(iF) CFMP_feature_2p_wells(iF) CFMP_feature_98p_wells(iF)...
        CFMP_feature_98p_wells(iF) CFMP_feature_2p_wells(iF)],[0.8 0.9 1],...
        'EdgeColor','none','linestyle','none');
    
    p2 = patch([numel(valid_control_ID)+10  nMovie nMovie numel(valid_control_ID)+10  numel(valid_control_ID)+10 ],...
        [CFMP_feature_2p_wells(iF) CFMP_feature_2p_wells(iF) CFMP_feature_98p_wells(iF)...
        CFMP_feature_98p_wells(iF) CFMP_feature_2p_wells(iF)],[0.8 0.9 1],...
        'EdgeColor','none','linestyle','none');
    
    hold on;
    line2 = plot([1 numel(valid_control_ID)],[CFMP_feature_2p_wells(iF) CFMP_feature_2p_wells(iF)],...
        '-.','color',[0.5 0.5 1]/1.02,'linewidth',2);
    line3 = plot([1 numel(valid_control_ID)],[CFMP_feature_50p_wells(iF) CFMP_feature_50p_wells(iF)],...
        '-','color',[0.5 0.5 1]/1.2,'linewidth',2);
    line4 = plot([1 numel(valid_control_ID)],[CFMP_feature_98p_wells(iF) CFMP_feature_98p_wells(iF)],...
        ':','color',[0.5 0.5 1]/1.02,'linewidth',2);
    
    
    
    line2 = plot([numel(valid_control_ID)+10 nMovie],[CFMP_feature_2p_wells(iF) CFMP_feature_2p_wells(iF)],...
        '-.','color',[0.5 0.5 1]/1.02,'linewidth',2);
    line3 = plot([numel(valid_control_ID)+10 nMovie],[CFMP_feature_50p_wells(iF) CFMP_feature_50p_wells(iF)],...
        '-','color',[0.5 0.5 1]/1.2,'linewidth',2);
    line4 = plot([numel(valid_control_ID)+10 nMovie],[CFMP_feature_98p_wells(iF) CFMP_feature_98p_wells(iF)],...
        ':','color',[0.5 0.5 1]/1.02,'linewidth',2);
   
    legend([line2,line3,line4],'2% of Control','Median of Control','98% of Control','FontSize',14);
    title(textcell{iF},'fontsize',15);
    set(gca, 'FontSize',14);
     
    % for samples
    % plot median for each one
    
    
     if iF==2    
   Control_thisFeature = squeeze(Control_ChMP_feature_thisML_ch1(iF,6,valid_control_ID));
    Positive_thisFeature = squeeze(Positive_ChMP_feature_thisML_ch1(iF,6,:));
     else
   Control_thisFeature = squeeze(Control_ChMP_feature_thisML_ch1(iF,5,valid_control_ID));
    Positive_thisFeature = squeeze(Positive_ChMP_feature_thisML_ch1(iF,5,:));
      end
    
    
    Y = Control_thisFeature;
    X = 1 :numel(Y);
    Y = Y';
    Y = Y(:);
    X_in = X(Y>= CFMP_feature_2p_wells(iF) & Y<= CFMP_feature_98p_wells(iF));
    Y_in = Y(Y>= CFMP_feature_2p_wells(iF) & Y<= CFMP_feature_98p_wells(iF));
    X_out = X(Y< CFMP_feature_2p_wells(iF) | Y> CFMP_feature_98p_wells(iF));
    Y_out = Y(Y< CFMP_feature_2p_wells(iF) | Y> CFMP_feature_98p_wells(iF));
        
    Y1=Y;
    try
    plot(X_in,Y_in,'.','MarkerSize',10);
    end
    try
    plot(X_out, Y_out,'*','MarkerSize',10);
    end
    
    axis([0 numel(valid_control_ID)+10+30 ...
        min(CFMP_feature_50p_wells(iF)+1.1*(min(Y)-CFMP_feature_50p_wells(iF)),CFMP_feature_50p_wells(iF)+2.5*(CFMP_feature_5p_wells(iF)-CFMP_feature_50p_wells(iF)))...
        max(CFMP_feature_50p_wells(iF)+1.1*(max(Y)-CFMP_feature_50p_wells(iF)),CFMP_feature_50p_wells(iF)+2.5*(CFMP_feature_95p_wells(iF)-CFMP_feature_50p_wells(iF)))]);
    
    for iX = 1:numel(X_out)
        iXX = X_out(iX);
        ID_str = ML_ID_cell{2}{valid_control_ID(iXX)};
        ID_str(1) = 'c';
        ID_str = [ID_str(1:2) '-' ID_str(3:end)];
        control_hits_count = control_hits_count+1;
        control_hits_pool(control_hits_count) = valid_control_ID(iXX);

        if(Y_out(iX)>CFMP_feature_50p_wells(iF))
        text(iXX-0.8, Y_out(iX)+(max(Y)-min(Y))*0.03, ID_str,'color','k','FontSize',14);
        else
        text(iXX-0.8, Y_out(iX)-(max(Y)-min(Y))*0.03, ID_str,'color','k','FontSize',14);
        end
    end
    
    
    if iF==2    
   Positive_thisFeature = squeeze(Positive_ChMP_feature_thisML_ch1(iF,6,:));
      else
   Positive_thisFeature = squeeze(Positive_ChMP_feature_thisML_ch1(iF,5,:));
       end
    
    Y = Positive_thisFeature;
    X = numel(valid_control_ID)+10+ (1:numel(Y));
    Y = Y';
    Y = Y(:);
    Y2 = Y;
    X_in = X(Y>= CFMP_feature_2p_wells(iF) & Y<= CFMP_feature_98p_wells(iF));
    Y_in = Y(Y>= CFMP_feature_2p_wells(iF) & Y<= CFMP_feature_98p_wells(iF));
    X_out = X(Y< CFMP_feature_2p_wells(iF) | Y> CFMP_feature_98p_wells(iF));
    Y_out = Y(Y< CFMP_feature_2p_wells(iF) | Y> CFMP_feature_98p_wells(iF));
        
    
    plot(X_in,Y_in,'r.','MarkerSize',10);
    plot(X_out, Y_out,'r*','MarkerSize',10);
    
    axis([0 numel(valid_control_ID)+10+30 ...
        min(CFMP_feature_50p_wells(iF)+1.1*(min([Y1; Y2])-CFMP_feature_50p_wells(iF)),CFMP_feature_50p_wells(iF)+2.5*(CFMP_feature_5p_wells(iF)-CFMP_feature_50p_wells(iF)))...
        max(CFMP_feature_50p_wells(iF)+1.1*(max([Y1; Y2])-CFMP_feature_50p_wells(iF)),CFMP_feature_50p_wells(iF)+2.5*(CFMP_feature_95p_wells(iF)-CFMP_feature_50p_wells(iF)))]);
     
    
     HitsList=[];
    for iX = 1:numel(X_out)
        iXX = X_out(iX) - (numel(valid_control_ID)+10);
        hits_count = hits_count+1;
        all_hits_pool(hits_count) = iXX;

        ID_str = ML_ID_cell{1}{iXX};
        ID_str(1) = 'p';
        ID_str = [ID_str(1:2) '-' ID_str(3:end)];
        HitsList = [HitsList,',',ID_str];
        if(Y_out(iX)>CFMP_feature_50p_wells(iF))
            text(X_out(iX)-0.8, Y_out(iX)+(max(Y)-min(Y))*0.03, ID_str,'color','m','FontSize',14);
        else
            text(X_out(iX)-0.8, Y_out(iX)-(max(Y)-min(Y))*0.03, ID_str,'color','m','FontSize',14);
        end
    end
    xlabel([num2str(numel(X_out)),' hits: ', HitsList(2:end)]);

    saveas(h1,[Group_ROOT_DIR,filesep,'Feature_',num2str(iF),'_Control_Positive.jpg']);
    saveas(h1,[Group_ROOT_DIR,filesep,'Feature_',num2str(iF),'_Control_Positive.tiff']);
    saveas(h1,[Group_ROOT_DIR,filesep,'Feature_',num2str(iF),'_Control_Positive.fig']);
    
    % for samples
    % plot median for each one
    
    
    h5=figure(5);
     
    if( iF==2    )
        Control_all = Control_ChMP_feature_thisML_ch1(iF,6,:);
    else
        Control_all = Control_ChMP_feature_thisML_ch1(iF,5,:);
    end
    
    Control_all(find(Control_Nuclueus <= 60 | Control_Nuclueus> 360))=nan;
    Control_all = Control_all(1:192);
    Control_all_matrix = nan(24,8);
    Control_all_matrix(:) = Control_all(:);
    
    plot_color_dot_from_matrix(flipud((Control_all_matrix)'),25);
    
    set(gca,'XTick',0.5:24-0.5);
    for iW = 1 : 24
      X_Tick_cell{iW} = num2str(iW);
    end
    set(gca,'XTickLabel',X_Tick_cell,'FontSize',14);

    set(gca,'YTick',0.5 :8-0.5);
     
    for iM = 1 :8
      Y_Tick_cell{iM} = ['c1-' char('A'+ 8 - iM)];
    end
   
    
    set(gca,'YTickLabel',Y_Tick_cell,'FontSize',14);
    CCCCC=  colormap;
    plot([0 0],[8.1 round(8*2)],'w','linewidth',10);
    plot([24.1 24*2],[0 0],'w','linewidth',10);
   
    saveas(h5,[Group_ROOT_DIR,filesep,'FeatureDots_',num2str(iF),'_control.jpg']);
    saveas(h5,[Group_ROOT_DIR,filesep,'FeatureDots_',num2str(iF),'_control.tiff']);
    saveas(h5,[Group_ROOT_DIR,filesep,'FeatureDots_',num2str(iF),'_control.fig']);
    
    
    h6=figure(6);
    if(iF==2)
        Positive_all = Positive_ChMP_feature_thisML_ch1(iF,6,:);
    else
        Positive_all = Positive_ChMP_feature_thisML_ch1(iF,5,:);
    end
    
    
    Positive_all_matrix = nan(5,6);
    Positive_all_matrix(:) = Positive_all(:);
    
    plot_color_dot_from_matrix(flipud((Positive_all_matrix)'),25);
    
    set(gca,'XTick',0.5:5-0.5);
    for iW = 1 : 5
      X_Tick_cell{iW} = num2str(iW+19);
    end
    set(gca,'XTickLabel',X_Tick_cell);

    set(gca,'YTick',0.5 :6-0.5);
     
    Y_Tick_cell = {'p3-P','p3-O','p2-P','p2-O','p1-P','p1-O'};
    
    set(gca,'YTickLabel',Y_Tick_cell);
    
    plot([0 0],[6.1 round(8*2)],'w','linewidth',10);
    plot([5.1 5*2],[0 0],'w','linewidth',10);
    colormap(CCCCC);    
   
    saveas(h6,[Group_ROOT_DIR,filesep,'FeatureDots_',num2str(iF),'_positive.jpg']);
    saveas(h6,[Group_ROOT_DIR,filesep,'FeatureDots_',num2str(iF),'_positive.tiff']);
    saveas(h6,[Group_ROOT_DIR,filesep,'FeatureDots_',num2str(iF),'_positive.fig']);
  
    

end