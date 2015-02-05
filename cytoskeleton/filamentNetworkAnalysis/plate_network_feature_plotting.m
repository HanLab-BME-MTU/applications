
% %% % plot the this group result
% h3 = figure(3);
% plot(Group_Pool_branch_duration_array, Group_Pool_branch_vif_mean_intensity,'.');
% xlabel(['Branch Duration,sample size:',num2str(numel(Group_Pool_branch_duration_array))],'Fontsize',13);
% ylabel('Branch Vim mean Int','Fontsize',13);
% title('Branch: Duration vs Vim Mean Int','Fontsize',13);
% set(gca,'fontsize',13)
% saveas(h3,[Group_ROOT_DIR,filesep,'EachBranch_Duration_vs_Vim.fig']);
% saveas(h3,[Group_ROOT_DIR,filesep,'EachBranch_Duration_vs_Vim.tif']);
% print(h3,'-depsc',[Group_ROOT_DIR,filesep,'EachBranch_Duration_vs_Vim.eps']);
% 
% 
% size(Group_Pool_protrusion_vif_mean_intensity)
% 
% h4 = figure(4);
% [h_p,bin] = hist(Group_Pool_protrusion_vif_mean_intensity,0.0:0.03:1.1);
% [h_r,bin] = hist(Group_Pool_retraction_vif_mean_intensity,0.0:0.03:1.1);
% bar([reshape(bin, 1, numel(bin)); reshape(bin, 1, numel(bin))]',[reshape(h_p, 1, numel(h_p));reshape(h_r, 1, numel(h_r))]');
% legend('Protrusion Region Vim','Retraction Region Vim');
% axis([0.00 0.30 0 max([h_p h_r])])
% set(gca,'fontsize',13)
% saveas(h4,[Group_ROOT_DIR,filesep,'Protrusion_vs_Retraction.fig']);
% saveas(h4,[Group_ROOT_DIR,filesep,'Protrusion_vs_Retraction.tif']);
% print(h4,'-depsc',[Group_ROOT_DIR,filesep,'Protrusion_vs_Retraction.eps']);
% 
% % mean vim vs speed
% h7 = figure(7); hold off;
% %  plot(Group_Pool_branch_number_mean, Group_Pool_Travel_Speed,'o');
% for i = 1 : length(Group_Pool_whole_cell_vif_mean_intensity)
%     % text(Group_Pool_branch_number_mean(i)-0.02, Group_Pool_Travel_Speed(i)+0.2, num2str(i), 'color',colorarray(i,:))
%     % hold on;plot(Group_Pool_branch_number_mean(i), Group_Pool_Travel_Speed(i), 'o', 'color',colorarray(i,:))
%     hold on;plot(Group_Pool_branch_number_mean(i), Group_Pool_Travel_Speed(i), 'bo','linewidth',2,'markersize',7);
%     
% end
% xlabel(['Average Branch Number, Sample Size:',num2str(numel(Group_Pool_branch_number_mean))],'FontSize',13);
% ylabel('Cell Speed','FontSize',13);
% title(['Branchness vs Speed, ',...
%     ' Correlation: ',...
%     num2str(corr(Group_Pool_branch_number_mean', ...
%     Group_Pool_Travel_Speed'),'%1.2f')],'FontSize',13);
% % axis([0 10 0 24]);
% set(gca,'fontsize',13);
% saveas(h7,[Group_ROOT_DIR,filesep,'Branchness_vs_Speed.fig']);
% saveas(h7,[Group_ROOT_DIR,filesep,'Branchness_vs_Speed.tif']);
% print(h7,'-depsc',[Group_ROOT_DIR,filesep,'Branchness_vs_Speed.eps']);
% 
% 
% h8 = figure(8);hold off;
% % plot(Group_Pool_branch_number_mean, Group_Pool_whole_cell_vif_mean_intensity,'.');
% for i = 1 : length(Group_Pool_whole_cell_vif_mean_intensity)
%     % text(Group_Pool_branch_number_mean(i)-0.02, Group_Pool_whole_cell_vif_mean_intensity(i)+10, num2str(i), 'color',colorarray(i,:));
%     hold on;
%     %     plot(Group_Pool_branch_number_mean(i), Group_Pool_whole_cell_vif_mean_intensity(i), '.', 'color',colorarray(i,:),'MarkerSize',16)
%     plot(Group_Pool_branch_number_mean(i), Group_Pool_whole_cell_vif_mean_intensity(i), 'bo','linewidth',2,'markersize',7);
% end
% xlabel(['Branch Mean Number, Sample Size:',num2str(numel(Group_Pool_branch_number_mean))],'Fontsize',13);
% ylabel('Cell Vim mean Int','Fontsize',13);
% title(['Branch number vs Vim, correlation: ',...
%     num2str(corr(Group_Pool_branch_number_mean', ...
%     Group_Pool_whole_cell_vif_mean_intensity'),'%1.2f') ],'Fontsize',13);
% set(gca,'fontsize',13);
% saveas(h8,[Group_ROOT_DIR,filesep,'Branchness_vs_Vim.fig']);
% saveas(h8,[Group_ROOT_DIR,filesep,'Branchness_vs_Vim.tif']);
% print(h8,'-depsc',[Group_ROOT_DIR,filesep,'Branchness_vs_Vim.eps']);
% 
% 
% h8_axis = axis;
% 
% h18 = figure(18);hold off;
% % plot(Group_Pool_branch_number_mean, Group_Pool_whole_cell_vif_mean_intensity,'.');
% for i = 1 : length(Group_Pool_whole_cell_vif_mean_intensity)
%     %     text(Group_Pool_thresholded_branch_number_mean(i)-0.02, Group_Pool_whole_cell_vif_mean_intensity(i)+10, num2str(i), 'color',colorarray(i,:));
%     hold on;
%     %     plot(Group_Pool_thresholded_branch_number_mean(i), Group_Pool_whole_cell_vif_mean_intensity(i), '.', 'color',colorarray(i,:),'MarkerSize',16);
%     plot(Group_Pool_thresholded_branch_number_mean(i), Group_Pool_whole_cell_vif_mean_intensity(i), 'bo','linewidth',2,'markersize',7);
% end
% corr_b_v_v = corrcoef(Group_Pool_whole_cell_vif_mean_intensity',Group_Pool_thresholded_branch_number_mean');
% display(['Corr for Branchness vs Vim: ',num2str(corr_b_v_v(1,2))]);
% 
% xlabel(['Branch Mean Number, Sample Size:',num2str(numel(Group_Pool_thresholded_branch_number_mean))],'Fontsize',13);
% ylabel('Cell Vim mean Int','Fontsize',13);
% title({['Branch number vs Vim, for branches larger than ',num2str(T_branchsize),' pixels'],['duration longer than ',num2str(T_branchduration),' frames'],...
%     ['Correlation:',num2str(corr_b_v_v(1,2))]},'Fontsize',13);
% axis(h8_axis);
% set(gca,'fontsize',13);
% saveas(h18,[Group_ROOT_DIR,filesep,'Branchness_vs_Vim_Thresholded.fig']);
% saveas(h18,[Group_ROOT_DIR,filesep,'Branchness_vs_Vim_Thresholded.tif']);
% print(h18,'-depsc',[Group_ROOT_DIR,filesep,'Branchness_vs_Vim_Thresholded.eps']);
% 
% 
% 
% %%
% try
%     %%
%     h26 = figure(26);hold off;
%     % plot(Group_Pool_fila_branch_orientation_pool_std, Group_Pool_whole_cell_vif_mean_intensity,'.');
%     for i = 1 : length(Group_Pool_whole_cell_vif_mean_intensity)
%         % text(Group_Pool_fila_branch_orientation_pool_std(i)-0.02, Group_Pool_whole_cell_vif_mean_intensity(i)+10, num2str(i), 'color',colorarray(i,:));
%         hold on;
%         % plot(Group_Pool_fila_branch_orientation_pool_std(i), Group_Pool_whole_cell_vif_mean_intensity(i), '.', 'color',colorarray(i,:),'MarkerSize',16);
%         plot(Group_Pool_fila_branch_orientation_pool_std(i), Group_Pool_whole_cell_vif_mean_intensity(i), 'bo','linewidth',2,'markersize',7);
%         
%     end
%     % axis([0.55, 1.0, 350 950]);
%     xlabel('Filament Orientation along Branch, Standard Deviation','Fontsize',13);
%     ylabel('Cell Vim mean Int','Fontsize',13);
%     title({'Filament Orientation Scatterness vs Vim Mean Int',...
%         ['Correlation: ',...
%         num2str(corr(Group_Pool_fila_branch_orientation_pool_std(~isnan(Group_Pool_fila_branch_orientation_pool_std)&~isnan(Group_Pool_whole_cell_vif_mean_intensity))', ...
%         Group_Pool_whole_cell_vif_mean_intensity(~isnan(Group_Pool_fila_branch_orientation_pool_std)&~isnan(Group_Pool_whole_cell_vif_mean_intensity))'),'%1.2f'), ...
%         ', Sample Size:',num2str(numel(Group_Pool_fila_branch_orientation_pool_std(~isnan(Group_Pool_fila_branch_orientation_pool_std)&~isnan(Group_Pool_whole_cell_vif_mean_intensity))))]},'Fontsize',13);
%     set(gca,'fontsize',13);
%     saveas(h26,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_branch_vs_Vim.fig']);
%     saveas(h26,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_branch_vs_Vim.tif']);
%     
%     print(h26,'-depsc',[Group_ROOT_DIR,filesep,'FilaOrient_wrt_branch_vs_Vim.eps']);
%     
%     
%     %%
%     
%     h36 = figure(36);hold off;
%     % plot(Group_Pool_fila_branch_orientation_pool_std, Group_Pool_whole_cell_vif_mean_intensity,'.');
%     for i = 1 : length(Group_Pool_branch_number_mean)
%         % text(Group_Pool_fila_branch_orientation_pool_std(i)-0.02, Group_Pool_branch_number_mean(i)+0.2, num2str(i), 'color',colorarray(i,:));
%         hold on;
%         %         plot(Group_Pool_fila_branch_orientation_pool_std(i), Group_Pool_branch_number_mean(i), '.', 'color',colorarray(i,:),'MarkerSize',16);
%         plot(Group_Pool_fila_branch_orientation_pool_std(i), Group_Pool_branch_number_mean(i), 'bo','linewidth',2,'markersize',7);
%     end
%     xlabel('Filament Orientation along Branch, Standard Deviation','Fontsize',13);
%     ylabel('Branch Mean Number ','Fontsize',13);
%     title({'Filament Orient along Branch Scatterness vs Branch Mean Num',...
%         ['Correlation: ',...
%         num2str(corr(Group_Pool_fila_branch_orientation_pool_std(~isnan(Group_Pool_fila_branch_orientation_pool_std)&~isnan(Group_Pool_branch_number_mean))', ...
%         Group_Pool_branch_number_mean(~isnan(Group_Pool_fila_branch_orientation_pool_std)&~isnan(Group_Pool_branch_number_mean))'),'%1.2f'),...
%         ', Sample Size:',num2str(numel(Group_Pool_fila_branch_orientation_pool_std(~isnan(Group_Pool_fila_branch_orientation_pool_std)&~isnan(Group_Pool_branch_number_mean))))]},'Fontsize',13);
%     axis([0.6 0.9 0 10]);
%     set(gca,'fontsize',13);
%     
%     saveas(h36,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_branch_vs_Branchness.fig']);
%     saveas(h36,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_branch_vs_Branchness.tif']);
%     print(h36,'-depsc',[Group_ROOT_DIR,filesep,'FilaOrient_wrt_branch_vs_Branchness.eps']);
%     
%     
% end
% 
% try
%     h46 = figure(46);hold off;
%     [h_result, bin_result] = hist(Group_Pool_fila_branch_orientation_pool,-pi/2+pi/36:pi/18:pi/2+pi/36);
%     bar(h_result/sum(h_result));
%     
%     set(gca,...
%         'xlim',[0.5 19.5],...
%         'xtick',[0.5:9:19.5],...
%         'xticklabel',{'-p/2'  '0' 'p/2' },...
%         'fontname','symbol',...
%         'fontsize',13)
%     
%     title('Filament Orientation along Branch Orientaion Distribution','Fontsize',13,'fontname','Helvetica');
%     set(gca,'fontsize',13);
%     saveas(h46,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_branchorient_distribution.fig']);
%     saveas(h46,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_branchorient_distribution.tif']);
%     
%     print(h46,'-depsc',[Group_ROOT_DIR,filesep,'FilaOrient_wrt_branchorient_distribution.eps']);
%     
%     
%     
%     h56 = figure(56);hold off;
%     [h_result, bin_result] = hist(Group_Pool_fila_branch_orientation_pool_fast,-pi/2+pi/36:pi/18:pi/2+pi/36);
%     bar(h_result/sum(h_result));
%     
%     set(gca,...
%         'xlim',[0.5 19.5],...
%         'xtick',[0.5:9:19.5],...
%         'ylim',[0.0 0.15],...
%         'ytick',[0.0:0.05:0.15],...
%         'yticklabel',{'0.0'  '0.05' '0.10' '0.15'},...
%         'xticklabel',{'-p/2'  '0' 'p/2' },...
%         'fontname','symbol',...
%         'fontsize',13)
%     title('Filament Orientation along Branch Orientaion Distribution(fast cell)','Fontsize',13,'fontname','Helvetica');
%     set(gca,'fontsize',13);
%     saveas(h56,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_branchorient_distribution_fastcell.fig']);
%     saveas(h56,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_branchorient_distribution_fastcell.tif']);
%     print(h56,'-depsc',[Group_ROOT_DIR,filesep,'FilaOrient_wrt_branchorient_distribution_fastcell.eps']);
%     
%     
%     h66 = figure(66);hold off;
%     [h_result, bin_result] = hist(Group_Pool_fila_branch_orientation_pool_slow,-pi/2+pi/36:pi/18:pi/2+pi/36);
%     bar(h_result/sum(h_result));
%     
%     set(gca,...
%         'xlim',[0.5 19.5],...
%         'xtick',[0.5:9:19.5],...
%         'ylim',[0.0 0.15],...
%         'ytick',[0.0:0.05:0.15],...
%         'yticklabel',{'0.0'  '0.05' '0.10' '0.15'},...
%         'xticklabel',{'-p/2'  '0' 'p/2' },...
%         'fontname','symbol',...
%         'fontsize',13)
%     set(gca,'fontsize',13);
%     title('Filament Orientation along Branch Orientaion Distribution(slow cell)','Fontsize',13,'fontname','Helvetica');
%     
%     saveas(h66,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_branchorient_distribution_slowcell.fig']);
%     saveas(h66,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_branchorient_distribution_slowcell.tif']);
%     print(h66,'-depsc',[Group_ROOT_DIR,filesep,'FilaOrient_wrt_branchorient_distribution_slowcell.eps']);
%     
%     
%     %%
%     
%     %%
%     h27 = figure(27);hold off;
%     % plot(Group_Pool_fila_trajectory_orientation_pool_std, Group_Pool_whole_cell_vif_mean_intensity,'.');
%     for i = 1 : length(Group_Pool_whole_cell_vif_mean_intensity)
%         % text(Group_Pool_fila_trajectory_orientation_pool_std(i)-0.02, Group_Pool_whole_cell_vif_mean_intensity(i)+10, num2str(i), 'color',colorarray(i,:));
%         hold on;
%         % plot(Group_Pool_fila_trajectory_orientation_pool_std(i), Group_Pool_whole_cell_vif_mean_intensity(i), '.', 'color',colorarray(i,:),'MarkerSize',16);
%         plot(Group_Pool_fila_trajectory_orientation_pool_std(i), Group_Pool_whole_cell_vif_mean_intensity(i), 'bo','linewidth',2,'markersize',7);
%         
%     end
%     % axis([0.55, 1.0, 350 950]);
%     xlabel('Filament Orientation along Path, Standard Deviation','Fontsize',13);
%     ylabel('Cell Vim mean Int','Fontsize',13);
%     title({'Filament Orientation along Path Scatterness vs Vim Mean Int',...
%         ['Correlation: ',...
%         num2str(corr(Group_Pool_fila_trajectory_orientation_pool_std(~isnan(Group_Pool_fila_trajectory_orientation_pool_std)&~isnan(Group_Pool_whole_cell_vif_mean_intensity))', ...
%         Group_Pool_whole_cell_vif_mean_intensity(~isnan(Group_Pool_fila_trajectory_orientation_pool_std)&~isnan(Group_Pool_whole_cell_vif_mean_intensity))'),'%1.2f'),...
%         ', Sample Size:',num2str(numel(Group_Pool_fila_branch_orientation_pool_std(~isnan(Group_Pool_fila_trajectory_orientation_pool_std)&~isnan(Group_Pool_whole_cell_vif_mean_intensity))))]},'Fontsize',13);
%     set(gca,'fontsize',13);
%     saveas(h27,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_path_vs_Vim.fig']);
%     saveas(h27,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_path_vs_Vim.tif']);
%     print(h27,'-depsc',[Group_ROOT_DIR,filesep,'FilaOrient_wrt_path_vs_Vim.eps']);
%     
%     
%     %%
% end
% 
% try
%     h37 = figure(37);hold off;
%     % plot(Group_Pool_fila_trajectory_orientation_pool_std, Group_Pool_whole_cell_vif_mean_intensity,'.');
%     for i = 1 : length(Group_Pool_branch_number_mean)
%         % text(Group_Pool_fila_trajectory_orientation_pool_std(i)-0.02, Group_Pool_branch_number_mean(i)+0.2, num2str(i), 'color',colorarray(i,:));
%         hold on;
%         %         plot(Group_Pool_fila_trajectory_orientation_pool_std(i), Group_Pool_branch_number_mean(i), '.', 'color',colorarray(i,:),'MarkerSize',16);
%         plot(Group_Pool_fila_trajectory_orientation_pool_std(i), Group_Pool_branch_number_mean(i),...
%             'bo','linewidth',2,'markersize',7);
%     end
%     xlabel('Filament Orientation along Cell Movement, Standard Deviation','Fontsize',13);
%     ylabel('Branch Mean Number ','Fontsize',13);
%     title({'Filament Orient along Cell Path Scatterness vs Branch Mean Num',...
%         ['Correlation: ',...
%         num2str(corr(Group_Pool_fila_trajectory_orientation_pool_std(~isnan(Group_Pool_fila_trajectory_orientation_pool_std)&~isnan(Group_Pool_branch_number_mean))', ...
%         Group_Pool_branch_number_mean(~isnan(Group_Pool_fila_trajectory_orientation_pool_std)&~isnan(Group_Pool_branch_number_mean))'),'%1.2f') ,...
%         ', Sample Size:',num2str(numel(Group_Pool_fila_trajectory_orientation_pool_std(~isnan(Group_Pool_fila_trajectory_orientation_pool_std)&~isnan(Group_Pool_branch_number_mean))))]},'Fontsize',13);
%     set(gca,'fontsize',13);
%     saveas(h37,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_path_vs_Branchness.fig']);
%     saveas(h37,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_path_vs_Branchness.tif']);
%     
%     print(h37,'-depsc',[Group_ROOT_DIR,filesep,'FilaOrient_wrt_path_vs_Branchness.eps']);
%     
%     
% end
% 
% try
%     h57 = figure(57);hold off;
%     [h_result, bin_result] = hist(Group_Pool_fila_trajectory_orientation_pool,-pi/2+pi/36:pi/18:pi/2+pi/36);
%     bar(h_result/sum(h_result));
%     
%     set(gca,...
%         'xlim',[0.5 19.5],...
%         'xtick',[0.5:9:19.5],...
%         'ylim',[0.0 0.15],...
%         'ytick',[0.0:0.05:0.15],...
%         'yticklabel',{'0.0'  '0.05' '0.10' '0.15'},...
%         'xticklabel',{'-p/2'  '0' 'p/2' },...
%         'fontname','symbol',...
%         'fontsize',13)
%     
%     title('Filament Orientation along Cell Movement Distribution','Fontsize',13,'fontname','Helvetica');
%     
%     saveas(h57,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_path_distribution.fig']);
%     saveas(h57,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_path_distribution.tif']);
%     
%     print(h57,'-depsc',[Group_ROOT_DIR,filesep,'FilaOrient_wrt_path_distribution.eps']);
%     
%     h57 = figure(57);hold off;
%     [h_result, bin_result] = hist(Group_Pool_fila_trajectory_orientation_pool_fast,-pi/2+pi/36:pi/18:pi/2+pi/36);
%     bar(h_result/sum(h_result));
%     
%     set(gca,...
%         'xlim',[0.5 19.5],...
%         'xtick',[0.5:9:19.5],...
%         'ylim',[0.0 0.15],...
%         'ytick',[0.0:0.05:0.15],...
%         'yticklabel',{'0.0'  '0.05' '0.10' '0.15'},...
%         'xticklabel',{'-p/2'  '0' 'p/2' },...
%         'fontname','symbol',...
%         'fontsize',13)
%     title('Filament Orientation along Cell Movement Distribution(fast cell)','Fontsize',13,'fontname','Helvetica');
%     
%     saveas(h57,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_path_distribution_fast.fig']);
%     saveas(h57,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_path_distribution_fast.tif']);
%     
%     print(h57,'-depsc',[Group_ROOT_DIR,filesep,'FilaOrient_wrt_path_distribution_fast.eps']);
%     
%     
%     h67 = figure(67);hold off;
%     [h_result, bin_result] = hist(Group_Pool_fila_trajectory_orientation_pool_slow,-pi/2+pi/36:pi/18:pi/2+pi/36);
%     bar(h_result/sum(h_result));
%     
%     set(gca,...
%         'xlim',[0.5 19.5],...
%         'xtick',[0.5:9:19.5],...
%         'ylim',[0.0 0.15],...
%         'ytick',[0.0:0.05:0.15],...
%         'yticklabel',{'0.0'  '0.05' '0.10' '0.15'},...
%         'xticklabel',{'-p/2'  '0' 'p/2' },...
%         'fontname','symbol',...
%         'fontsize',13)
%     title('Filament Orientation along Cell Movement Distribution(slow cells)','Fontsize',13,'fontname','Helvetica');
%     set(gca,'fontsize',13);
%     saveas(h67,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_path_distribution_slow.fig']);
%     saveas(h67,[Group_ROOT_DIR,filesep,'FilaOrient_wrt_path_distribution_slow.tif']);
%     print(h67,'-depsc',[Group_ROOT_DIR,filesep,'FilaOrient_wrt_path_distribution_slow.eps']);
%     
% end
% %%
% h28 = figure(28);hold off;
% % plot(Group_Pool_branch_cellmovement_std, Group_Pool_whole_cell_vif_mean_intensity,'.');
% for i = 1 : length(Group_Pool_whole_cell_vif_mean_intensity)
%     % text(Group_Pool_branch_cellmovement_std(i)-0.02, Group_Pool_whole_cell_vif_mean_intensity(i)+10, num2str(i), 'color',colorarray(i,:));
%     hold on;
%     % plot(Group_Pool_branch_cellmovement_std(i), Group_Pool_whole_cell_vif_mean_intensity(i), '.', 'color',colorarray(i,:),'MarkerSize',16);
%     plot(Group_Pool_branch_cellmovement_std(i), Group_Pool_whole_cell_vif_mean_intensity(i), 'bo','linewidth',2,'markersize',7);
%     
% end
% % axis([0.55, 1.0, 350 950]);
% xlabel('Branch Orientation along Cell Movement, Standard Deviation','Fontsize',13);
% ylabel('Cell Vim mean Int','Fontsize',13);
% title({'Branch Orientation Scatterness vs Vim',...
%     ['Correlation: ',...
%     num2str(corr(Group_Pool_branch_cellmovement_std(~isnan(Group_Pool_branch_cellmovement_std)&~isnan(Group_Pool_whole_cell_vif_mean_intensity))', ...
%     Group_Pool_whole_cell_vif_mean_intensity(~isnan(Group_Pool_branch_cellmovement_std)&~isnan(Group_Pool_whole_cell_vif_mean_intensity))'),'%1.2f') ,...
%     ', Sample Size:',num2str(numel(Group_Pool_whole_cell_vif_mean_intensity(~isnan(Group_Pool_whole_cell_vif_mean_intensity)&~isnan(Group_Pool_whole_cell_vif_mean_intensity))))]},'Fontsize',13);
% set(gca,'fontsize',13);
% saveas(h28,[Group_ROOT_DIR,filesep,'BranchOrient_vs_Vim.fig']);
% saveas(h28,[Group_ROOT_DIR,filesep,'BranchOrient_vs_Vim.tif']);
% print(h28,'-depsc',[Group_ROOT_DIR,filesep,'BranchOrient_vs_Vim.eps']);
% 
% 
% %%
% 
% h38 = figure(38);hold off;
% % plot(Group_Pool_branch_cellmovement_std, Group_Pool_whole_cell_vif_mean_intensity,'.');
% for i = 1 : length(Group_Pool_branch_number_mean)
%     % text(Group_Pool_branch_cellmovement_std(i)-0.02, Group_Pool_branch_number_mean(i)+0.2, num2str(i), 'color',colorarray(i,:));
%     hold on;
%     %     plot(Group_Pool_branch_cellmovement_std(i), Group_Pool_branch_number_mean(i), '.', 'color',colorarray(i,:),'MarkerSize',16);
%     plot(Group_Pool_branch_cellmovement_std(i), Group_Pool_branch_number_mean(i), 'bo','linewidth',2,'markersize',7);
% end
% xlabel('Branch Orientation along Cell Movement, Standard Deviation','Fontsize',13);
% ylabel('Branch Mean Number ','Fontsize',13);
% title({'Branch Orientation Scatterness vs Branch Mean Number',...
%     ['Correlation: ',...
%     num2str(corr(Group_Pool_branch_cellmovement_std(~isnan(Group_Pool_branch_cellmovement_std)&~isnan(Group_Pool_branch_number_mean))', ...
%     Group_Pool_branch_number_mean(~isnan(Group_Pool_branch_cellmovement_std)&~isnan(Group_Pool_branch_number_mean))'),'%1.2f') ,...
%     ', Sample Size:',num2str(numel(Group_Pool_branch_number_mean(~isnan(Group_Pool_branch_number_mean)&~isnan(Group_Pool_branch_cellmovement_std))))]},'Fontsize',13);
% set(gca,'fontsize',13);
% saveas(h38,[Group_ROOT_DIR,filesep,'BranchOrient_vs_Branchness.fig']);
% saveas(h38,[Group_ROOT_DIR,filesep,'BranchOrient_vs_Branchness.tif']);
% print(h38,'-depsc',[Group_ROOT_DIR,filesep,'BranchOrient_vs_Branchness.eps']);
% 
% h48 = figure(48);hold off;
% % plot(Group_Pool_branch_cellmovement_std, Group_Pool_whole_cell_vif_mean_intensity,'.');
% for i = 1 : numel(Group_Pool_whole_cell_vif_mean_intensity)
%     text(Group_Pool_Travel_Speed(i)-0.02, ...
%         Group_Pool_whole_cell_vif_mean_intensity(i)+0.02, ...
%         Identifier_cell{i}, 'color',colorarray(i,:));
%     hold on;
%     %     plot(Group_Pool_Travel_Speed(i), Group_Pool_whole_cell_vif_mean_intensity(i), '.', 'color',colorarray(i,:),'MarkerSize',16);
%     plot(Group_Pool_Travel_Speed(i), Group_Pool_whole_cell_vif_mean_intensity(i), 'bo','linewidth',2,'markersize',7);
% end
% xlabel('Cell Speed','Fontsize',13);
% ylabel('Vim Level ','Fontsize',13);
% title({'Cell Speed vs Vim Level',...
%     ['Correlation: ',...
%     num2str(corr(Group_Pool_Travel_Speed', ...
%     Group_Pool_whole_cell_vif_mean_intensity'),'%1.2f'),...
%     ', Sample Size:',num2str(numel(Group_Pool_Travel_Speed))]},'Fontsize',13);
% set(gca,'fontsize',13);
% saveas(h48,[Group_ROOT_DIR,filesep,'Speed_vs_Vim_color.fig']);
% saveas(h48,[Group_ROOT_DIR,filesep,'Speed_vs_Vim_color.tif']);
% print(h48,'-depsc',[Group_ROOT_DIR,filesep,'Speed_vs_Vim_color.eps']);
% 
% 
% 
% h58 = figure(58);hold off;
% 
% % plot(Group_Pool_branch_cellmovement_std, Group_Pool_whole_cell_vif_mean_intensity,'.');
% for i = 1 : length(Group_Pool_whole_cell_vif_mean_intensity)
%     % text(Group_Pool_Travel_Speed(i)-0.02, Group_Pool_whole_cell_vif_mean_intensity(i)+0.2, num2str(i), 'color',colorarray(i,:));
%     if(Group_Pool_Cell_Marked_Frame_Number(i)>15 && Group_Pool_CompletedFrame_last(i)>60)
%         hold on;
%         %     plot(Group_Pool_Travel_Speed(i), Group_Pool_whole_cell_vif_mean_intensity(i), '.', 'color',colorarray(i,:),'MarkerSize',16);
%         plot(Group_Pool_Travel_Speed(i), Group_Pool_whole_cell_vif_mean_intensity(i), 'bo','linewidth',2,'markersize',7);
%     end
% end
% xlabel('Cell Speed','Fontsize',13);
% ylabel('Vim Level ','Fontsize',13);
% title({'Cell Speed vs Vim Level (cell with more frames, with late hours)',...
%     ['Correlation: ',...
%     num2str(corr(Group_Pool_Travel_Speed(Group_Pool_Cell_Marked_Frame_Number>15 & Group_Pool_CompletedFrame_last>60)', ...
%     Group_Pool_whole_cell_vif_mean_intensity(Group_Pool_Cell_Marked_Frame_Number>15 & Group_Pool_CompletedFrame_last>60)'),'%1.2f'),...
%     ', Sample Size:',num2str(numel(Group_Pool_Travel_Speed(Group_Pool_Cell_Marked_Frame_Number>15 & Group_Pool_CompletedFrame_last>60)))]},'Fontsize',13);
% set(gca,'fontsize',13);
% saveas(h58,[Group_ROOT_DIR,filesep,'Speed_vs_Vim_moreframes.fig']);
% saveas(h58,[Group_ROOT_DIR,filesep,'Speed_vs_Vim_moreframes.tif']);
% print(h58,'-depsc',[Group_ROOT_DIR,filesep,'Speed_vs_Vim_moreframes.eps']);
% 
% 
% 
% %%
% %
% h9 = figure(9);
% plot(Group_Pool_branch_number_mean_pat, Group_Pool_branch_vif_mean_intensity,'.');
% xlabel('Branch Mean Number','fontsize',13);
% ylabel('Branch Vim mean Int','fontsize',13);
% title('Branch number vs Vim','fontsize',13);
% 
% set(gca,'fontsize',13);
% saveas(h9,[Group_ROOT_DIR,filesep,'Branchness_vs_BranchVim.fig']);
% saveas(h9,[Group_ROOT_DIR,filesep,'Branchness_vs_BranchVim.tif']);
% print(h9,'-depsc',[Group_ROOT_DIR,filesep,'Branchness_vs_BranchVim.eps']);
% 
% 
% h10 = figure(10);hold off;
% % plot(Group_Pool_branch_number_mean, Group_Pool_whole_cell_vif_mean_intensity,'.');
% for i = 1 : length(Group_Pool_whole_cell_vim_totalamount_mean)
%     %     text(Group_Pool_branch_number_mean(i)-0.02, Group_Pool_whole_cell_vim_totalamount_mean(i)+10, num2str(i), 'color',colorarray(i,:));
%     hold on;
%     %     plot(Group_Pool_branch_number_mean(i), Group_Pool_whole_cell_vim_totalamount_mean(i), '.', 'color',colorarray(i,:))
%     plot(Group_Pool_branch_number_mean(i), Group_Pool_whole_cell_vim_totalamount_mean(i), 'bo','linewidth',2,'markersize',7);
%     
% end
% xlabel('Branch Mean Number','fontsize',13);
% ylabel('Cell Vim Total Int','fontsize',13);
% title({['Branch number vs Vim Total(Sample number: ', num2str(length(Group_Pool_whole_cell_vim_totalamount_mean)),')'],...
%     ['Correlation: ',num2str(corr(Group_Pool_branch_number_mean',Group_Pool_whole_cell_vim_totalamount_mean'))]},'fontsize',13);
% set(gca,'fontsize',13);
% 
% saveas(h10,[Group_ROOT_DIR,filesep,'Branchness_vs_VimTotal.fig']);
% saveas(h10,[Group_ROOT_DIR,filesep,'Branchness_vs_VimTotal.tif']);
% print(h10,'-depsc',[Group_ROOT_DIR,filesep,'Branchness_vs_VimTotal.eps']);
% 
% 
% % h20 = figure(20);hold off;
% % % plot(Group_Pool_branch_number_mean, Group_Pool_whole_cell_vif_mean_intensity,'.');
% % for i = 1 : length(Group_Pool_whole_cell_vif_mean_intensity)
% % % text(Group_Pool_branch_number_mean(i)-0.02, Group_Pool_whole_cell_vif_mean_intensity(i)+10, num2str(i), 'color',colorarray(i,:));
% % hold on;
% % plot(Group_Pool_branch_number_mean(i), Group_Pool_whole_cell_vif_mean_intensity(i), '.', 'color',colorarray(i,:),'MarkerSize',15)
% % end
% % xlabel('Branch Mean Number','Fontsize',13);
% % ylabel('Cell Vim Mean Int','Fontsize',13);
% % title(['Branch number vs Vim Mean,  correlation: ',...
% %     num2str(corr(Group_Pool_branch_number_mean', ...
% %     Group_Pool_whole_cell_vif_mean_intensity'),'%1.2f') ],'Fontsize',13);
% % saveas(h20,[Group_ROOT_DIR,filesep,'Branchness_vs_VimMean.fig']);
% % saveas(h20,[Group_ROOT_DIR,filesep,'Branchness_vs_VimMean.tif']);
% %
% % h30 = figure(30);hold off;
% % % plot(Group_Pool_branch_number_mean, Group_Pool_whole_cell_vif_mean_intensity,'.');
% % for i = 1 : length(Group_Pool_whole_cell_vif_mean_intensity)
% % % text(Group_Pool_branch_number_mean(i)-0.02, Group_Pool_whole_cell_vif_mean_intensity(i)+10, num2str(i), 'color',colorarray(i,:));
% % hold on;
% % plot(Group_Pool_Travel_Speed(i), Group_Pool_whole_cell_vif_mean_intensity(i), '.', 'color',colorarray(i,:),'MarkerSize',15)
% % end
% % xlabel('Cell Speed','Fontsize',13);
% % ylabel('Cell Vim Mean Int','Fontsize',13);
% % title(['Cell Speed vs Vim Mean,  correlation: ',...
% %     num2str(corr(Group_Pool_Travel_Speed', ...
% %     Group_Pool_whole_cell_vif_mean_intensity'),'%1.2f') ],'Fontsize',13);
% % saveas(h30,[Group_ROOT_DIR,filesep,'Speed_vs_VimMean.fig']);
% % saveas(h30,[Group_ROOT_DIR,filesep,'Speed_vs_VimMean.tif']);
% 
% 
% h11 = figure(11);hold off;
% % plot(Group_Pool_branch_number_mean, Group_Pool_whole_cell_vif_mean_intensity,'.');
% for i = 1 : length(Group_Pool_tracked_branch_d_frames)
%     %     text(Group_Pool_tracked_branch_d_frames(i)-0.02, Group_Pool_whole_cell_vim_totalamount_mean(i)+10, num2str(i), 'color',colorarray(i,:));
%     hold on;
%     %     plot(Group_Pool_tracked_branch_d_frames(i), Group_Pool_whole_cell_vim_totalamount_mean(i), '.', 'color',colorarray(i,:))
%     plot(Group_Pool_tracked_branch_d_frames(i), Group_Pool_whole_cell_vim_totalamount_mean(i), 'bo','linewidth',2,'markersize',7);
% end
% xlabel('Tracked Branchs per Frame','fontsize',13);
% ylabel('Cell Vim Total Int','fontsize',13);
% title({['Tracked Branches per Frame vs Vim Total(Sample number: ', num2str(length(Group_Pool_whole_cell_vim_totalamount_mean)),')'],...
%     ['Correlation: ',num2str(corr(Group_Pool_tracked_branch_d_frames',Group_Pool_whole_cell_vim_totalamount_mean'))]},'fontsize',13);
% set(gca,'fontsize',13);
% 
% saveas(h11,[Group_ROOT_DIR,filesep,'TrackedBranchness_vs_VimTotal.fig']);
% saveas(h11,[Group_ROOT_DIR,filesep,'TrackedBranchness_vs_VimTotal.tif']);
% print(h11,'-depsc',[Group_ROOT_DIR,filesep,'TrackedBranchness_vs_VimTotal.eps']);
% 
% 
% h12 = figure(12);
% plot(Group_Pool_whole_cell_vif_mean_intensity_pat, Group_Pool_branch_duration_array, '.');
% ylabel('Branch Duration','fontsize',13);
% xlabel('Cell Vim mean Int','fontsize',13);
% title('Cell Vim vs Branch Duration','fontsize',13);
% set(gca,'fontsize',13);
% 
% saveas(h12,[Group_ROOT_DIR,filesep,'BranchDuration_vs_CellVim.fig']);
% saveas(h12,[Group_ROOT_DIR,filesep,'BranchDuration_vs_CellVim.tif']);
% print(h12,'-depsc',[Group_ROOT_DIR,filesep,'BranchDuration_vs_CellVim.eps']);
% 
% 
% h13 = figure(13);
% plot(Group_Pool_branch_duration_array(find(Group_Pool_branch_size_mean>T_branchsize & Group_Pool_branch_duration_array>T_branchduration)), ...
%     Group_Pool_branch_vif_mean_intensity(find(Group_Pool_branch_size_mean>T_branchsize & Group_Pool_branch_duration_array>T_branchduration)),'.');
% xlabel('Branch Duration','fontsize',13);
% ylabel('Branch Vim mean Int','fontsize',13);
% title({['Branch: Duration vs Vim, for branches larger than ',num2str(T_branchsize),' pixels,'],['duration longer than ',num2str(T_branchduration),' frames']},'fontsize',13);
% set(gca,'fontsize',13);
% saveas(h13,[Group_ROOT_DIR,filesep,'EachBranch_Duration_vs_Vim_with_Thres.fig']);
% saveas(h13,[Group_ROOT_DIR,filesep,'EachBranch_Duration_vs_Vim_with_Thres.tif']);
% print(h13,'-depsc',[Group_ROOT_DIR,filesep,'EachBranch_Duration_vs_Vim_with_Thres.eps']);
% 
% save([Group_ROOT_DIR,filesep,'branch_analysis_group_results.mat']);
% 
% 
% %%
% branch_analysis_plotting_VimNms;
% branch_analysis_plotting_VimTotalNms;
% branch_analysis_plotting_VimFilaDen;
% branch_analysis_plotting_VimFilaTotal;
% 
% 
