
h108 = figure(108);close;
h118 = figure(118);close;
h118 = figure(128);close;
h148 = figure(148);close;
h158 = figure(158);close;

Group_Pool_whole_cell_vim_seg_mean = Group_Pool_whole_cell_vim_seg_mean(:)';

h108 = figure(108);hold off;
% plot(Group_Pool_branch_number_mean, Group_Pool_whole_cell_vif_mean_intensity,'.');
for i = 1 : length(Group_Pool_whole_cell_vif_mean_intensity)
    % text(Group_Pool_branch_number_mean(i)-0.02, Group_Pool_whole_cell_vif_mean_intensity(i)+10, num2str(i), 'color',colorarray(i,:));
    hold on;
    %     plot(Group_Pool_branch_number_mean(i), Group_Pool_whole_cell_vif_mean_intensity(i), '.', 'color',colorarray(i,:),'MarkerSize',16)
    plot(Group_Pool_branch_number_mean(i), Group_Pool_whole_cell_vim_seg_mean(i), 'bo','linewidth',2,'markersize',7);
end
xlabel(['Branch Mean Number, Sample Size:',num2str(numel(Group_Pool_branch_number_mean))],'Fontsize',13);
ylabel('Cell Vim mean filament density','Fontsize',13);
title(['Branch number vs Vim mean filament density, correlation: ',...
    num2str(corr(Group_Pool_branch_number_mean(~isnan(Group_Pool_branch_number_mean)&~isnan(Group_Pool_whole_cell_vim_seg_mean))', ...
    Group_Pool_whole_cell_vim_seg_mean(~isnan(Group_Pool_branch_number_mean)&~isnan(Group_Pool_whole_cell_vim_seg_mean))'),'%1.2f') ],'Fontsize',13);
set(gca,'fontsize',13);
saveas(h108,[Group_ROOT_DIR,filesep,'Branchness_vs_VimFilaDen.fig']);
saveas(h108,[Group_ROOT_DIR,filesep,'Branchness_vs_VimFilaDen.tif']);
print(h108,'-depsc',[Group_ROOT_DIR,filesep,'Branchness_vs_VimFilaDen.eps']);


h108_axis = axis;

h118 = figure(118);hold off;
% plot(Group_Pool_branch_number_mean, Group_Pool_whole_cell_vif_mean_intensity,'.');
for i = 1 : length(Group_Pool_whole_cell_vif_mean_intensity)
    %     text(Group_Pool_thresholded_branch_number_mean(i)-0.02, Group_Pool_whole_cell_vif_mean_intensity(i)+10, num2str(i), 'color',colorarray(i,:));
    hold on;
    %     plot(Group_Pool_thresholded_branch_number_mean(i), Group_Pool_whole_cell_vif_mean_intensity(i), '.', 'color',colorarray(i,:),'MarkerSize',16);
    plot(Group_Pool_thresholded_branch_number_mean(i), Group_Pool_whole_cell_vim_seg_mean(i), 'bo','linewidth',2,'markersize',7);
end
corr_b_v_v = corrcoef(Group_Pool_whole_cell_vim_seg_mean(~isnan(Group_Pool_branch_number_mean)&~isnan(Group_Pool_whole_cell_vim_seg_mean))',...
    Group_Pool_thresholded_branch_number_mean(~isnan(Group_Pool_branch_number_mean)&~isnan(Group_Pool_whole_cell_vim_seg_mean))');
display(['Corr for Branchness vs Vim mean filament density: ',num2str(corr_b_v_v(1,2))]);

xlabel(['Branch Mean Number, Sample Size:',num2str(numel(Group_Pool_thresholded_branch_number_mean))],'Fontsize',13);
ylabel('Cell Vim mean filament density','Fontsize',13);
title({['Branch number vs Vim mean filament density, for branches larger than ',num2str(T_branchsize),' pixels'],['duration longer than ',num2str(T_branchduration),' frames'],...
    ['Correlation:',num2str(corr_b_v_v(1,2))]},'Fontsize',13);
axis(h108_axis);
set(gca,'fontsize',13);
saveas(h118,[Group_ROOT_DIR,filesep,'Branchness_vs_VimFilaDen_Thresholded.fig']);
saveas(h118,[Group_ROOT_DIR,filesep,'Branchness_vs_VimFilaDen_Thresholded.tif']);
print(h118,'-depsc',[Group_ROOT_DIR,filesep,'Branchness_vs_VimFilaDen_Thresholded.eps']);

%%
h128 = figure(128);hold off;
% plot(Group_Pool_branch_cellmovement_std, Group_Pool_whole_cell_vif_mean_intensity,'.');
for i = 1 : length(Group_Pool_whole_cell_vim_seg_mean)
    % text(Group_Pool_branch_cellmovement_std(i)-0.02, Group_Pool_whole_cell_vif_mean_intensity(i)+10, num2str(i), 'color',colorarray(i,:));
    hold on;
    % plot(Group_Pool_branch_cellmovement_std(i), Group_Pool_whole_cell_vif_mean_intensity(i), '.', 'color',colorarray(i,:),'MarkerSize',16);
    plot(Group_Pool_branch_cellmovement_std(i), Group_Pool_whole_cell_vim_seg_mean(i), 'bo','linewidth',2,'markersize',7);
    
end
% axis([0.55, 1.0, 350 950]);
xlabel('Branch Orientation along Cell Movement, Standard Deviation','Fontsize',13);
ylabel('Cell Vim mean filament density','Fontsize',13);
title({'Branch Orientation Scatterness vs Vim mean filament density',...
    ['Correlation: ',...
    num2str(corr(Group_Pool_branch_cellmovement_std(~isnan(Group_Pool_branch_cellmovement_std)&~isnan(Group_Pool_whole_cell_vim_seg_mean))', ...
    Group_Pool_whole_cell_vim_seg_mean(~isnan(Group_Pool_branch_cellmovement_std)&~isnan(Group_Pool_whole_cell_vim_seg_mean))'),'%1.2f') ,...
    ', Sample Size:',num2str(numel(Group_Pool_whole_cell_vif_mean_intensity(~isnan(Group_Pool_branch_cellmovement_std)&~isnan(Group_Pool_whole_cell_vim_seg_mean))))]},'Fontsize',13);
set(gca,'fontsize',13);
saveas(h128,[Group_ROOT_DIR,filesep,'BranchOrient_vs_VimFilaDen.fig']);
saveas(h128,[Group_ROOT_DIR,filesep,'BranchOrient_vs_VimFilaDen.tif']);
print(h128,'-depsc',[Group_ROOT_DIR,filesep,'BranchOrient_vs_VimFilaDen.eps']);

%%
h148 = figure(148);hold off;
% plot(Group_Pool_branch_cellmovement_std, Group_Pool_whole_cell_vif_mean_intensity,'.');
for i = 1 : numel(Group_Pool_whole_cell_vif_mean_intensity)
    text(Group_Pool_Travel_Speed(i)-0.02, ...
        Group_Pool_whole_cell_vim_seg_mean(i)+0.002, ...
        Identifier_cell{i}, 'color',colorarray(i,:));
    hold on;
    %     plot(Group_Pool_Travel_Speed(i), Group_Pool_whole_cell_vif_mean_intensity(i), '.', 'color',colorarray(i,:),'MarkerSize',16);
    plot(Group_Pool_Travel_Speed(i), Group_Pool_whole_cell_vim_seg_mean(i), 'bo','linewidth',2,'markersize',7);
end
xlabel('Cell Speed','Fontsize',13);
ylabel('Vim Filament Density','Fontsize',13);
title({'Cell Speed vs Vim mean filament density',...
    ['Correlation: ',...
    num2str(corr(Group_Pool_Travel_Speed(~isnan(Group_Pool_Travel_Speed)&~isnan(Group_Pool_whole_cell_vim_seg_mean))', ...
    Group_Pool_whole_cell_vim_seg_mean(~isnan(Group_Pool_Travel_Speed)&~isnan(Group_Pool_whole_cell_vim_seg_mean))'),'%1.2f'),...
    ', Sample Size:',num2str(numel(Group_Pool_Travel_Speed))]},'Fontsize',13);
set(gca,'fontsize',13);
saveas(h148,[Group_ROOT_DIR,filesep,'Speed_vs_VimFilaDen_color.fig']);
saveas(h148,[Group_ROOT_DIR,filesep,'Speed_vs_VimFilaDen_color.tif']);
print(h148,'-depsc',[Group_ROOT_DIR,filesep,'Speed_vs_VimFilaDen_color.eps']);

h158 = figure(158);hold off;

% plot(Group_Pool_branch_cellmovement_std, Group_Pool_whole_cell_vif_mean_intensity,'.');
for i = 1 : length(Group_Pool_whole_cell_vim_seg_mean)
    % text(Group_Pool_Travel_Speed(i)-0.02, Group_Pool_whole_cell_vif_mean_intensity(i)+0.2, num2str(i), 'color',colorarray(i,:));
    if(Group_Pool_Cell_Marked_Frame_Number(i)>15 && Group_Pool_CompletedFrame_last(i)>60)
    hold on;
    %     plot(Group_Pool_Travel_Speed(i), Group_Pool_whole_cell_vif_mean_intensity(i), '.', 'color',colorarray(i,:),'MarkerSize',16);
    plot(Group_Pool_Travel_Speed(i), Group_Pool_whole_cell_vim_seg_mean(i), 'bo','linewidth',2,'markersize',7);
    end
end
xlabel('Cell Speed','Fontsize',13);
ylabel('Vim Filament Density ','Fontsize',13);

X = Group_Pool_Travel_Speed(Group_Pool_Cell_Marked_Frame_Number>15 & Group_Pool_CompletedFrame_last>60);
Y = Group_Pool_whole_cell_vim_seg_mean(Group_Pool_Cell_Marked_Frame_Number>15 & Group_Pool_CompletedFrame_last>60);

title({'Cell Speed vs Vim mean filament density (cell with more frames, with late hours)',...
    ['Correlation: ',...
    num2str(corr(X(~isnan(X)&~isnan(Y))', ...
    Y(~isnan(X)&~isnan(Y))'),'%1.2f'),...
    ', Sample Size:',num2str(numel(Group_Pool_Travel_Speed(Group_Pool_Cell_Marked_Frame_Number>15 & Group_Pool_CompletedFrame_last>60)))]},'Fontsize',13);
set(gca,'fontsize',13);
saveas(h158,[Group_ROOT_DIR,filesep,'Speed_vs_VimFilaDen_moreframes.fig']);
saveas(h158,[Group_ROOT_DIR,filesep,'Speed_vs_VimFilaDen_moreframes.tif']);
print(h158,'-depsc',[Group_ROOT_DIR,filesep,'Speed_vs_VimFilaDen_moreframes.eps']);
