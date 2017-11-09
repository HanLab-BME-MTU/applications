function plot_analysis_hist_from_twoMDs(MD_filename1, MD_filename2)
% function to plot the filament intensity and scale for two movies

% Input: MD_filename1: the file name of the movieData.mat file for movie 1
%        MD_filename2: the file name of the movieData.mat file for movie 2

% Output: None

% Created 07 2014 by Liya Ding, Matlab R2012b


try
    load(MD_filename1);
    MD1 = MD;
catch
    display('Invalid MD filename 1');
    return;
end

try
    load(MD_filename2);
    MD2 = MD;
catch
    display('Invalid MD filename 2');
    return;
end

if(numel(MD1.channels_)~=2 ||numel(MD2.channels_)~=2)
     display('Invalid channel number');   
    return;
end


%% pool the analysis together

feature_pool_1 = pool_filament_network_analysis_MD(MD1);

feature_pool_2 = pool_filament_network_analysis_MD(MD2);


%% scale plotting
[hist11,b11] = hist(feature_pool_1{1}.scale_per_filament_pool_all,1:0.05:2);
hist11 = hist11./(sum(hist11))*100;

[hist21,b21] = hist(feature_pool_2{1}.scale_per_filament_pool_all,1:0.05:2);
hist21 = hist21./(sum(hist21))*100;

[hist12,b12] = hist(feature_pool_1{2}.scale_per_filament_pool_all,1:0.05:2);
hist12 = hist12./(sum(hist12))*100;

[hist22,b22] = hist(feature_pool_2{2}.scale_per_filament_pool_all,1:0.05:2);
hist22 = hist22./(sum(hist22))*100;

h1 = figure(1);

subplot(221);bar(hist11);
axis([0.5 21.5 0 max([max(hist11),max(hist21)]) ]);
set(gca,'XTick',1:4:21);
set(gca,'XTickLabel',{'1.0','1.2','1.4','1.6','1.8','2.0'});

ylabel('Percent(%)')

title('MD1 Channel1 Scale');
subplot(222);bar(hist12);
axis([0.5 21.5 0 max([max(hist12),max(hist22)]) ]);
set(gca,'XTick',1:4:21);
set(gca,'XTickLabel',{'1.0','1.2','1.4','1.6','1.8','2.0'});
ylabel('Percent(%)')
title('MD1 Channel2 Scale');

subplot(223);bar(hist21);
axis([0.5 21.5 0 max([max(hist11),max(hist21)]) ]);
set(gca,'XTick',1:4:21);
set(gca,'XTickLabel',{'1.0','1.2','1.4','1.6','1.8','2.0'});

ylabel('Percent(%)')

title('MD2 Channel1 Scale');
subplot(224);bar(hist22);
axis([0.5 21.5 0 max([max(hist12),max(hist22)]) ]);
set(gca,'XTick',1:4:21);
set(gca,'XTickLabel',{'1.0','1.2','1.4','1.6','1.8','2.0'});
ylabel('Percent(%)')
title('MD2 Channel2 Scale');


%% mean intensity

mean_intensity_per_filament_pool_all1 = [...
    feature_pool_1{1}.mean_intensity_per_filament_pool_all;...
    feature_pool_2{1}.mean_intensity_per_filament_pool_all;];
 
   mean_intensity_per_filament_pool_all2 = [...
       feature_pool_1{2}.mean_intensity_per_filament_pool_all;...
       feature_pool_2{2}.mean_intensity_per_filament_pool_all;];

MN1 = min(mean_intensity_per_filament_pool_all1);
MX1 = max(mean_intensity_per_filament_pool_all1);

MN2 = min(mean_intensity_per_filament_pool_all2);
MX2 = max(mean_intensity_per_filament_pool_all2);

[hist11,b11] = hist(feature_pool_1{1}.mean_intensity_per_filament_pool_all,MN1:(MX1-MN1)/20:MX1);
hist11 = hist11./(sum(hist11))*100;
[hist21,b21] = hist(feature_pool_2{1}.mean_intensity_per_filament_pool_all,MN1:(MX1-MN1)/20:MX1);
hist21 = hist21./(sum(hist21))*100;

[hist12,b12] = hist(feature_pool_1{2}.mean_intensity_per_filament_pool_all,MN2:(MX2-MN2)/20:MX2);
hist12 = hist12./(sum(hist12))*100;
[hist22,b22] = hist(feature_pool_2{2}.mean_intensity_per_filament_pool_all,MN2:(MX2-MN2)/20:MX2);
hist22 = hist22./(sum(hist22))*100;

h2 = figure(2);
subplot(221);bar(hist11);
axis([0.5 21.5 0 max([max(hist11),max(hist21)]) ]);
set(gca,'XTick',1:4:21);
set(gca,'XTickLabel',{num2str(MN1,'%.0f'),...
    num2str(MN1+(MX1-MN1)/20*4,'%.0f'),...
    num2str(MN1+(MX1-MN1)/20*8,'%.0f'),...
    num2str(MN1+(MX1-MN1)/20*12,'%.0f'),...
    num2str(MN1+(MX1-MN1)/20*16,'%.0f'),...
    num2str(MN1+(MX1-MN1)/20*20,'%.0f')});

ylabel('Percent(%)')
title('MD1 Channel1 Filament Intenisty');

subplot(222);bar(hist12);
axis([0.5 21.5 0 max([max(hist12),max(hist22)]) ]);
set(gca,'XTick',1:4:21);
set(gca,'XTickLabel',{num2str(MN2,'%.0f'),...
    num2str(MN2+(MX2-MN2)/20*4,'%.0f'),...
    num2str(MN2+(MX2-MN2)/20*8,'%.0f'),...
    num2str(MN2+(MX2-MN2)/20*12,'%.0f'),...
    num2str(MN2+(MX2-MN2)/20*16,'%.0f'),...
    num2str(MN2+(MX2-MN2)/20*20,'%.0f')});

ylabel('Percent(%)')
title('MD1 Channel2 Filament Intenisty');

subplot(223);bar(hist21);
axis([0.5 21.5 0 max([max(hist11),max(hist21)]) ]);
set(gca,'XTick',1:4:21);
set(gca,'XTickLabel',{num2str(MN1,'%.0f'),...
    num2str(MN1+(MX1-MN1)/20*4,'%.0f'),...
    num2str(MN1+(MX1-MN1)/20*8,'%.0f'),...
    num2str(MN1+(MX1-MN1)/20*12,'%.0f'),...
    num2str(MN1+(MX1-MN1)/20*16,'%.0f'),...
    num2str(MN1+(MX1-MN1)/20*20,'%.0f')});
ylabel('Percent(%)')
title('MD2 Channel1 Filament Intenisty');

subplot(224);bar(hist22);
axis([0.5 21.5 0 max([max(hist12),max(hist22)]) ]);
set(gca,'XTick',1:4:21);
set(gca,'XTickLabel',{num2str(MN2,'%.0f'),...
    num2str(MN2+(MX2-MN2)/20*4,'%.0f'),...
    num2str(MN2+(MX2-MN2)/20*8,'%.0f'),...
    num2str(MN2+(MX2-MN2)/20*12,'%.0f'),...
    num2str(MN2+(MX2-MN2)/20*16,'%.0f'),...
    num2str(MN2+(MX2-MN2)/20*20,'%.0f')});


ylabel('Percent(%)')
title('MD2 Channel2 Filament Intenisty');


%%
%% mean fat filament intensity

mean_intensity_per_fat_filament_pool_all1 = [...
    feature_pool_1{1}.mean_intensity_per_fat_filament_pool_all;...
    feature_pool_2{1}.mean_intensity_per_fat_filament_pool_all;];
 
   mean_intensity_per_fat_filament_pool_all2 = [...
       feature_pool_1{2}.mean_intensity_per_fat_filament_pool_all;...
       feature_pool_2{2}.mean_intensity_per_fat_filament_pool_all;];

MN1 = min(mean_intensity_per_fat_filament_pool_all1);
MX1 = max(mean_intensity_per_fat_filament_pool_all1);

MN2 = min(mean_intensity_per_fat_filament_pool_all2);
MX2 = max(mean_intensity_per_fat_filament_pool_all2);

[hist11,b11] = hist(feature_pool_1{1}.mean_intensity_per_fat_filament_pool_all,MN1:(MX1-MN1)/20:MX1);
hist11 = hist11./(sum(hist11))*100;
[hist21,b21] = hist(feature_pool_2{1}.mean_intensity_per_fat_filament_pool_all,MN1:(MX1-MN1)/20:MX1);
hist21 = hist21./(sum(hist21))*100;

[hist12,b12] = hist(feature_pool_1{2}.mean_intensity_per_fat_filament_pool_all,MN2:(MX2-MN2)/20:MX2);
hist12 = hist12./(sum(hist12))*100;
[hist22,b22] = hist(feature_pool_2{2}.mean_intensity_per_fat_filament_pool_all,MN2:(MX2-MN2)/20:MX2);
hist22 = hist22./(sum(hist22))*100;

h3 = figure(3);
subplot(221);bar(hist11);
axis([0.5 21.5 0 max([max(hist11),max(hist21)]) ]);
set(gca,'XTick',1:4:21);
set(gca,'XTickLabel',{num2str(MN1,'%.0f'),...
    num2str(MN1+(MX1-MN1)/20*4,'%.0f'),...
    num2str(MN1+(MX1-MN1)/20*8,'%.0f'),...
    num2str(MN1+(MX1-MN1)/20*12,'%.0f'),...
    num2str(MN1+(MX1-MN1)/20*16,'%.0f'),...
    num2str(MN1+(MX1-MN1)/20*20,'%.0f')});

ylabel('Percent(%)')
title('MD1 Channel1 Fat Filament Int');

subplot(222);bar(hist12);
axis([0.5 21.5 0 max([max(hist12),max(hist22)]) ]);
set(gca,'XTick',1:4:21);
set(gca,'XTickLabel',{num2str(MN2,'%.0f'),...
    num2str(MN2+(MX2-MN2)/20*4,'%.0f'),...
    num2str(MN2+(MX2-MN2)/20*8,'%.0f'),...
    num2str(MN2+(MX2-MN2)/20*12,'%.0f'),...
    num2str(MN2+(MX2-MN2)/20*16,'%.0f'),...
    num2str(MN2+(MX2-MN2)/20*20,'%.0f')});

ylabel('Percent(%)')
title('MD1 Channel2 Fat Filament Int');

subplot(223);bar(hist21);
axis([0.5 21.5 0 max([max(hist11),max(hist21)]) ]);
set(gca,'XTick',1:4:21);
set(gca,'XTickLabel',{num2str(MN1,'%.0f'),...
    num2str(MN1+(MX1-MN1)/20*4,'%.0f'),...
    num2str(MN1+(MX1-MN1)/20*8,'%.0f'),...
    num2str(MN1+(MX1-MN1)/20*12,'%.0f'),...
    num2str(MN1+(MX1-MN1)/20*16,'%.0f'),...
    num2str(MN1+(MX1-MN1)/20*20,'%.0f')});
ylabel('Percent(%)')
title('MD2 Channel1 Fat Filament Int');

subplot(224);bar(hist22);
axis([0.5 21.5 0 max([max(hist12),max(hist22)]) ]);
set(gca,'XTick',1:4:21);
set(gca,'XTickLabel',{num2str(MN2,'%.0f'),...
    num2str(MN2+(MX2-MN2)/20*4,'%.0f'),...
    num2str(MN2+(MX2-MN2)/20*8,'%.0f'),...
    num2str(MN2+(MX2-MN2)/20*12,'%.0f'),...
    num2str(MN2+(MX2-MN2)/20*16,'%.0f'),...
    num2str(MN2+(MX2-MN2)/20*20,'%.0f')});

ylabel('Percent(%)')
title('MD2 Channel2 Fat Filament Int');


%% mean steerable filtering responce

mean_steerable_per_filament_pool_all1 = [...
    feature_pool_1{1}.mean_steerable_per_filament_pool_all;...
    feature_pool_2{1}.mean_steerable_per_filament_pool_all;];
 
   mean_steerable_per_filament_pool_all2 = [...
       feature_pool_1{2}.mean_steerable_per_filament_pool_all;...
       feature_pool_2{2}.mean_steerable_per_filament_pool_all;];

MN1 = min(mean_steerable_per_filament_pool_all1);
MX1 = max(mean_steerable_per_filament_pool_all1);

MN2 = min(mean_steerable_per_filament_pool_all2);
MX2 = max(mean_steerable_per_filament_pool_all2);

[hist11,b11] = hist(feature_pool_1{1}.mean_steerable_per_filament_pool_all,MN1:(MX1-MN1)/20:MX1);
hist11 = hist11./(sum(hist11))*100;
[hist21,b21] = hist(feature_pool_2{1}.mean_steerable_per_filament_pool_all,MN1:(MX1-MN1)/20:MX1);
hist21 = hist21./(sum(hist21))*100;

[hist12,b12] = hist(feature_pool_1{2}.mean_steerable_per_filament_pool_all,MN2:(MX2-MN2)/20:MX2);
hist12 = hist12./(sum(hist12))*100;
[hist22,b22] = hist(feature_pool_2{2}.mean_steerable_per_filament_pool_all,MN2:(MX2-MN2)/20:MX2);
hist22 = hist22./(sum(hist22))*100;

h4 = figure(4);
subplot(221);bar(hist11);
axis([0.5 21.5 0 max([max(hist11),max(hist21)]) ]);
set(gca,'XTick',1:4:21);
set(gca,'XTickLabel',{num2str(MN1,'%.0f'),...
    num2str(MN1+(MX1-MN1)/20*4,'%.0f'),...
    num2str(MN1+(MX1-MN1)/20*8,'%.0f'),...
    num2str(MN1+(MX1-MN1)/20*12,'%.0f'),...
    num2str(MN1+(MX1-MN1)/20*16,'%.0f'),...
    num2str(MN1+(MX1-MN1)/20*20,'%.0f')});

ylabel('Percent(%)')
title('MD1 Channel1 Filament Steerable Fitering Responce');

subplot(222);bar(hist12);
axis([0.5 21.5 0 max([max(hist12),max(hist22)]) ]);
set(gca,'XTick',1:4:21);
set(gca,'XTickLabel',{num2str(MN2,'%.0f'),...
    num2str(MN2+(MX2-MN2)/20*4,'%.0f'),...
    num2str(MN2+(MX2-MN2)/20*8,'%.0f'),...
    num2str(MN2+(MX2-MN2)/20*12,'%.0f'),...
    num2str(MN2+(MX2-MN2)/20*16,'%.0f'),...
    num2str(MN2+(MX2-MN2)/20*20,'%.0f')});

ylabel('Percent(%)')
title('MD1 Channel2 Filament Steerable Fitering Responce');

subplot(223);bar(hist21);
axis([0.5 21.5 0 max([max(hist11),max(hist21)]) ]);
set(gca,'XTick',1:4:21);
set(gca,'XTickLabel',{num2str(MN1,'%.0f'),...
    num2str(MN1+(MX1-MN1)/20*4,'%.0f'),...
    num2str(MN1+(MX1-MN1)/20*8,'%.0f'),...
    num2str(MN1+(MX1-MN1)/20*12,'%.0f'),...
    num2str(MN1+(MX1-MN1)/20*16,'%.0f'),...
    num2str(MN1+(MX1-MN1)/20*20,'%.0f')});
ylabel('Percent(%)')
title('MD2 Channel1 Filament Steerable Fitering Responce');

subplot(224);bar(hist22);
axis([0.5 21.5 0 max([max(hist12),max(hist22)]) ]);
set(gca,'XTick',1:4:21);
set(gca,'XTickLabel',{num2str(MN2,'%.0f'),...
    num2str(MN2+(MX2-MN2)/20*4,'%.0f'),...
    num2str(MN2+(MX2-MN2)/20*8,'%.0f'),...
    num2str(MN2+(MX2-MN2)/20*12,'%.0f'),...
    num2str(MN2+(MX2-MN2)/20*16,'%.0f'),...
    num2str(MN2+(MX2-MN2)/20*20,'%.0f')});


ylabel('Percent(%)')
title('MD2 Channel2 Filament Steerable Fitering Responce');


% 
% AA = [min([C_mt.output_feature.mean_intensity_per_fat_filament_pool P_mt.output_feature.mean_intensity_per_fat_filament_pool])...
%     max([C_mt.output_feature.mean_intensity_per_fat_filament_pool P_mt.output_feature.mean_intensity_per_fat_filament_pool]) 0 250];
% BB = [min([C_vif.output_feature.mean_intensity_per_fat_filament_pool P_vif.output_feature.mean_intensity_per_fat_filament_pool])...
%     max([C_vif.output_feature.mean_intensity_per_fat_filament_pool P_vif.output_feature.mean_intensity_per_fat_filament_pool]) 0 120];
% 
% figure;subplot(221);hist(C_mt.output_feature.mean_intensity_per_fat_filament_pool,AA(1):(AA(2)-AA(1))/20:AA(2));
% axis(AA);
% title('Control MT Mean Int(Fat)');
% subplot(222);hist(C_vif.output_feature.mean_intensity_per_fat_filament_pool,BB(1):(BB(2)-BB(1))/20:BB(2));
% axis(BB);
% title('Control VIF Mean Int(Fat)');
% 
% subplot(223);hist(P_mt.output_feature.mean_intensity_per_fat_filament_pool,AA(1):(AA(2)-AA(1))/20:AA(2));
% axis(AA);
% title('Positive MT Mean Int(Fat)');
% subplot(224);hist(P_vif.output_feature.mean_intensity_per_fat_filament_pool,BB(1):(BB(2)-BB(1))/20:BB(2));
% axis(BB);
% title('Positive VIF Mean Int(Fat)');
% 
