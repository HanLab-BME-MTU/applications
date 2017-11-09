% function network_feature_ML_cell = load_ML_network_correlation_analysis(ML,radius,figure_flag)
% % function to do network analysis for a whole movielist
% % Liya Ding, June, 2014
% %
% % Input:
% %   ML:     The movieList object loaded before running this function
% 
% % the number of movies
% 
% % if(~exist('radius','var'))
% %     radius=20;
% % end


if(~exist('figure_flag','var'))
    figure_flag=0;
end


movieNumber =  length(ML.movieDataFile_);
network_feature_ML_cell= cell(1,1);

for iM  = 1 : movieNumber
    
    clearvars -except 'movieNumber' 'network_feature_ML_cell'...
        'iM' 'ML' 'figure_flag' 'radius'
    
    close all;
    
    % load this movie
    load(ML.movieDataFile_{iM});
    % the number of channels
    display('======================================================================');
    display(['iM:', num2str(iM)]);
    
    display_msg_flag = 0; % display warning or not
    package_process_ind_script;
    
    for iChannel = 1 : numel(MD.channels_)
        outdir = [MD.processes_{indexFilamentSegmentationProcess}.outFilePaths_{iChannel},filesep,'analysis_results'];
           Channel_FilesNames = MD.channels_(iChannel).getImageFileNames(1:MD.nFrames_);
    
           filename_short_strs = uncommon_str_takeout(Channel_FilesNames);

        
        for iFrame = 1 : MD.nFrames_
         try
            load([outdir,filesep,'network_feature_ch_',num2str(iChannel),'_',filename_short_strs{iFrame},'.mat'],...
                'output_feature');
            
            network_feature{iChannel,iFrame} = output_feature;
         catch
                 network_feature{iChannel,iFrame} = [];   
         end
        end
       
    end
   
    try        
        network_feature_ML_cell{1, iM} = network_feature;        
    catch        
        network_feature_ML_cell{1, iM} = load_MD_network_for_analysis(MD,[],radius,[ones(5,1); zeros(11,1)]);
    end
end

ML_ROOT_DIR = ML.outputDirectory_;
% save([ML_ROOT_DIR,'\movieList_netwrok_analysis_output','(radius',num2str(radius),').0727.mat'],...
%     'network_feature_ML_cell');


%% plot the density out

MT_density_all_pool = [];
VIF_density_all_pool = [];
scrable_MT_density_all_pool = [];
scrable_VIF_density_all_pool = [];
corr_array = [];
corr_array_S_VIF=[];
corr_array_S_MT=[];
corr_array_SS=[];

for iM  = 1 :movieNumber
    load(ML.movieDataFile_{iM});
    
    MT_density_pool = [];
    VIF_density_pool = [];
    scrable_MT_density_pool=[];
    scrable_VIF_density_pool = [];    
    
    this_MD_f = network_feature_ML_cell{1, iM};
    nFrame = size(this_MD_f,2);
    
    for iF = 1 : nFrame
        MT_density = this_MD_f{1,1}.density_filament;
        VIF_density = this_MD_f{2,1}.density_filament;
        
        scrable_MT_density = this_MD_f{1,1}.scrabled_density_filament;
        scrable_VIF_density = this_MD_f{2,1}.scrabled_density_filament;        
        
        MT_cell_seg = this_MD_f{1,1}.Cell_Mask;
        VIF_cell_seg = this_MD_f{2,1}.Cell_Mask;
        Cell_Mask = MT_cell_seg.*VIF_cell_seg;
        
        Cell_Mask = MD.processes_{2}.loadChannelOutput(1,1);
        
         Cell_Mask(1:radius+1,:)=0;
        Cell_Mask(:,1:radius+1)=0;
        Cell_Mask(end-radius:end,:)=0;
        Cell_Mask(:,end-radius:end)=0;

        MT_density_pool = [MT_density_pool MT_density(Cell_Mask>0)];
        VIF_density_pool = [VIF_density_pool VIF_density(Cell_Mask>0)];
        scrable_MT_density_pool = [scrable_MT_density_pool scrable_MT_density(Cell_Mask>0)];
        scrable_VIF_density_pool = [scrable_VIF_density_pool scrable_VIF_density(Cell_Mask>0)];
    end
    
    h11= figure(11);hold off;
    
    A = MT_density_pool(1:5:end);
B = VIF_density_pool(1:5:end);
C = corrcoef(A(~isnan(A)&~isnan(B)),B(~isnan(A)&~isnan(B)));

    plot(MT_density_pool(1:50:end),VIF_density_pool(1:50:end),'.');
    title(['Movie-',num2str(iM),', MT and VIF densities Cross-correlation:', ...
        num2str(C(1,2)),'(radius: ',num2str(radius),')']);
    saveas(h11, [ML_ROOT_DIR,filesep,'density_correlation_movie',num2str(iM),'(radius',num2str(radius),').0727.tif']);
    saveas(h11, [ML_ROOT_DIR,filesep,'density_correlation_movie',num2str(iM),'(radius',num2str(radius),').0727.fig']);
    corr_array(iM) = C(1,2);
    
           A = MT_density_pool(1:5:end);
B = scrable_VIF_density_pool(1:5:end);
C = corrcoef(A(~isnan(A)&~isnan(B)),B(~isnan(A)&~isnan(B)));
 
    h12= figure(12);hold off;
    plot(MT_density_pool(1:50:end),scrable_VIF_density_pool(1:50:end),'.');
    title(['Movie-',num2str(iM),', MT and scrabled VIF densities Cross-correlation:', ...
        num2str(C(1,2)),'(radius: ',num2str(radius),')']);
    saveas(h12, [ML_ROOT_DIR,filesep,'scrable_VIF_density_correlation_movie',num2str(iM),'(radius',num2str(radius),').0727.tif']);
    saveas(h12, [ML_ROOT_DIR,filesep,'scrable_VIF_density_correlation_movie',num2str(iM),'(radius',num2str(radius),').0727.fig']);
    corr_array_S_VIF(iM) = C(1,2);
    
     A = scrable_MT_density_pool(1:5:end);
B = VIF_density_pool(1:5:end);
C = corrcoef(A(~isnan(A)&~isnan(B)),B(~isnan(A)&~isnan(B)));

    h13= figure(13);hold off;
    plot(scrable_MT_density_pool(1:50:end),VIF_density_pool(1:50:end),'.');
    title(['Movie-',num2str(iM),', scrabled MT and VIF densities Cross-correlation:', ...
        num2str(C(1,2)),'(radius: ',num2str(radius),')']);
    saveas(h13, [ML_ROOT_DIR,filesep,'scrable_MT_density_correlation_movie',num2str(iM),'(radius',num2str(radius),').0727.tif']);
    saveas(h13, [ML_ROOT_DIR,filesep,'scrable_MT_density_correlation_movie',num2str(iM),'(radius',num2str(radius),').0727.fig']);
    corr_array_S_MT(iM) = C(1,2);
    
     A = scrable_MT_density_pool(1:5:end);
B = scrable_VIF_density_pool(1:5:end);
C = corrcoef(A(~isnan(A)&~isnan(B)),B(~isnan(A)&~isnan(B)));
 h14= figure(14);hold off;
    plot(scrable_MT_density_pool(1:50:end),scrable_VIF_density_pool(1:50:end),'.');
    title(['Movie-',num2str(iM),', scrabled MT and VIF densities Cross-correlation:', ...
        num2str(C(1,2)),'(radius: ',num2str(radius),')']);
    saveas(h14, [ML_ROOT_DIR,filesep,'scrable_MT_VIF_density_correlation_movie',num2str(iM),'(radius',num2str(radius),').0727.tif']);
    saveas(h14, [ML_ROOT_DIR,filesep,'scrable_MT_VIF_density_correlation_movie',num2str(iM),'(radius',num2str(radius),').0727.fig']);
    corr_array_SS(iM) = C(1,2);
        
    MT_density_all_pool = [MT_density_all_pool; MT_density_pool];
    VIF_density_all_pool = [VIF_density_all_pool; VIF_density_pool];
    
    scrable_MT_density_all_pool = [scrable_MT_density_all_pool; scrable_MT_density_pool];
    scrable_VIF_density_all_pool = [scrable_VIF_density_all_pool; scrable_VIF_density_pool];
    
end

h21= figure(21);
plot(MT_density_all_pool(1:100:end),VIF_density_all_pool(1:100:end),'.');
title(['MT and VIF densities Cross-correlation: ', ...
    num2str(corr(MT_density_all_pool(1:50:end),VIF_density_all_pool(1:50:end))),'(radius: ',num2str(radius),')']);
saveas(h21, [ML_ROOT_DIR,filesep,'scrable_VIF_density_correlation_movie',num2str(iM),'(radius',num2str(radius),').0727.tif']);
saveas(h21, [ML_ROOT_DIR,filesep,'scrable_VIF_density_correlation_movie',num2str(iM),'(radius',num2str(radius),').0727.fig']);
AA = axis;

h22= figure(22);hold off;
plot(MT_density_all_pool(1:100:end),scrable_VIF_density_all_pool(1:100:end),'.');
title(['Movie-',num2str(iM),', MT and scrabled VIF densities Cross-correlation:', ...
    num2str(corr(MT_density_all_pool(1:50:end),scrable_VIF_density_all_pool(1:50:end))),'(radius: ',num2str(radius),')']);
axis(AA);
saveas(h22, [ML_ROOT_DIR,filesep,'scrable_VIF_density_correlation_movie',num2str(iM),'(radius',num2str(radius),').0727.tif']);
saveas(h22, [ML_ROOT_DIR,filesep,'scrable_VIF_density_correlation_movie',num2str(iM),'(radius',num2str(radius),').0727.fig']);

h23= figure(23);hold off;
plot(scrable_MT_density_all_pool(1:100:end),VIF_density_all_pool(1:100:end),'.');
title(['Movie-',num2str(iM),', scrabled MT and VIF densities Cross-correlation:', ...
    num2str(corr(scrable_MT_density_all_pool(1:50:end),VIF_density_all_pool(1:50:end))),'(radius: ',num2str(radius),')']);
axis(AA);
saveas(h23, [ML_ROOT_DIR,filesep,'scrable_MT_density_correlation_movie',num2str(iM),'(radius',num2str(radius),').0727.tif']);
saveas(h23, [ML_ROOT_DIR,filesep,'scrable_MT_density_correlation_movie',num2str(iM),'(radius',num2str(radius),').0727.fig']);


h24= figure(24);hold off;
plot(scrable_MT_density_all_pool(1:100:end),scrable_VIF_density_all_pool(1:100:end),'.');
title(['Movie-',num2str(iM),', scrabled MT and VIF densities Cross-correlation:', ...
    num2str(corr(scrable_MT_density_all_pool(1:50:end),scrable_VIF_density_all_pool(1:50:end))),'(radius: ',num2str(radius),')']);
axis(AA);
saveas(h24, [ML_ROOT_DIR,filesep,'scrable_MT_VIF_density_correlation_movie',num2str(iM),'(radius',num2str(radius),').0727.tif']);
saveas(h24, [ML_ROOT_DIR,filesep,'scrable_MT_VIF_density_correlation_movie',num2str(iM),'(radius',num2str(radius),').0727.fig']);

corr_ML = [corr(MT_density_all_pool,VIF_density_all_pool);
    corr(MT_density_all_pool,scrable_VIF_density_all_pool);
    corr(scrable_MT_density_all_pool,VIF_density_all_pool);
    corr(scrable_MT_density_all_pool,scrable_VIF_density_all_pool);];

save([ML_ROOT_DIR,filesep,'ML_corr_radius_',num2str(radius),'.0727.mat'],...
    'corr_array','corr_array_S_VIF',...
    'corr_array_S_MT','corr_array_SS',...
    'corr_ML');
