function network_feature_ML_cell = load_ML_network_analysis(ML,radius,figure_flag, save_everything_flag,...
    feature_flag, vimscreen_flag, wholeML_flag)
% function to do network analysis for a whole movielist
% Liya Ding, June, 2014
%
% Input:
%   ML:     The movieList object loaded before running this function
%           radius: the definition of neighborhood
%           figure_flag: 1 to plot histgrams, 0 not to
%           save_everything_flag: 1 to save the plot, 0 not to
%           feature_flag: the flag for controlling which feature to
%                         calculate, a vector of 25 bits, 1 for need to
%                         calculate, 0 for skip it to save time.
%           vimscreen_flag: flag for vim screen
%           wholeML_flag: 1 if want to put everything in this ML together,
%                         by default 0


if(~exist('radius','var'))
    radius = 20;
end

if(~exist('figure_flag','var'))
    figure_flag = 0;
end

if(~exist('save_everything_flag','var'))
    save_everything_flag = 0;
end

if(~exist('feature_flag','var'))
    feature_flag = ones(25,1);
end

if(~exist('vimscreen_flag','var'))
    vimscreen_flag = 0;
end

if(~exist('wholeML_flag','var'))
    wholeML_flag = 0;
end


feature_index = feature_index(:);  
if(numel(feature_index)<25)
    feature_index = [feature_index; zeros(25-numel(feature_index), 1)];
end

%%
movieNumber =  length(ML.movieDataFile_);
network_feature_ML_cell = cell(1,1);
network_feature_ML_wholepool_cell = cell(1,1);

for iM  = 1 :movieNumber
    
    clearvars -except 'movieNumber' 'network_feature_ML_cell'...
        'network_feature_ML_wholepool_cell'...
        'iM' 'ML' 'figure_flag' 'radius'...
        'save_everything_flag'...
        'vimscreen_flag' 'wholeML_flag'

    close all;
    
    % load this movie
    MD = MovieData.load(ML.movieDataFile_{iM});
    % the number of channels
    display('======================================================================');
    display(['iM:', num2str(iM)]);
    
    % if want to keep everything together into a cell, get the output
    if(wholeML_flag==1)
       [network_feature_MD_cell,network_feature_MD_wholepool_cell] = load_MD_network_for_analysis...
        (MD,[],radius,figure_flag, save_everything_flag,feature_flag,vimscreen_flag);
    else
        % if just want to run every movie, run without bring back
        % output(everything is saved to harddisk already)
        load_MD_network_for_analysis...
        (MD,[],radius,figure_flag, save_everything_flag,feature_flag,vimscreen_flag);
    end
    
    if(wholeML_flag==1)
        network_feature_ML_cell{1, iM} = network_feature_MD_cell;
        network_feature_ML_wholepool_cell{1, iM} = network_feature_MD_wholepool_cell;
    end
end

if(wholeML_flag==1)
    ML_ROOT_DIR = ML.outputDirectory_;
    save([ML_ROOT_DIR, filesep,'movieList_netwrok_analysis_output','(radius',num2str(radius),').mat'],...
        'network_feature_ML_cell','network_feature_ML_wholepool_cell');
    
end
    
% %% plot the density out
% 
% MT_density_all_pool = [];
% VIF_density_all_pool = [];
% scrable_MT_density_all_pool = [];
% scrable_VIF_density_all_pool = [];
% 
% corr_array = [];
% corr_array_S_VIF=[];
% corr_array_S_MT=[];
% corr_array_SS=[];
% 
% for iM  = 1 : movieNumber
%     
%     MT_density_pool = [];
%     VIF_density_pool = [];
%     scrable_MT_density_pool=[];
%     scrable_VIF_density_pool = [];    
%     
%     this_MD_f = network_feature_ML_cell{1, iM};
%     nFrame = size(this_MD_f,2);
%     
%     for iF = 1 : nFrame
%         MT_density = this_MD_f{1,1}.density_filament;
%         VIF_density = this_MD_f{2,1}.density_filament;
%         
%         scrable_MT_density = this_MD_f{1,1}.scrabled_density_filament;
%         scrable_VIF_density = this_MD_f{2,1}.scrabled_density_filament;
%         
%         MT_cell_seg = this_MD_f{1,1}.Cell_Mask;
%         VIF_cell_seg = this_MD_f{2,1}.Cell_Mask;
%         Cell_Mask = MT_cell_seg.*VIF_cell_seg;
%         
%         MT_density_pool = [MT_density_pool MT_density(Cell_Mask>0)];
%         VIF_density_pool = [VIF_density_pool VIF_density(Cell_Mask>0)];
%         scrable_MT_density_pool = [scrable_MT_density_pool scrable_MT_density(Cell_Mask>0)];
%         scrable_VIF_density_pool = [scrable_VIF_density_pool scrable_VIF_density(Cell_Mask>0)];
%     end
%     
%     h11= figure(11);hold off;
%     plot(MT_density_pool(1:50:end),VIF_density_pool(1:50:end),'.');
%     title(['Movie-',num2str(iM),', MT and VIF densities Cross-correlation:', ...
%         num2str(corr(MT_density_pool,VIF_density_pool)),'(radius: ',num2str(radius),')']);
%     saveas(h11, [ML_ROOT_DIR,filesep,'density_correlation_movie',num2str(iM),'(radius',num2str(radius),').tif']);
%     saveas(h11, [ML_ROOT_DIR,filesep,'density_correlation_movie',num2str(iM),'(radius',num2str(radius),').fig']);
%     corr_array(iM) = corr(MT_density_pool,VIF_density_pool);
%     
%         
%     h12= figure(12);hold off;
%     plot(MT_density_pool(1:50:end),scrable_VIF_density_pool(1:50:end),'.');
%     title(['Movie-',num2str(iM),', MT and scrabled VIF densities Cross-correlation:', ...
%         num2str(corr(MT_density_pool,scrable_VIF_density_pool)),'(radius: ',num2str(radius),')']);
%     saveas(h12, [ML_ROOT_DIR,filesep,'scrable_VIF_density_correlation_movie',num2str(iM),'(radius',num2str(radius),').tif']);
%     corr_array_S_VIF(iM) = corr(MT_density_pool,scrable_VIF_density_pool);
%     
%     h13= figure(13);hold off;
%     plot(scrable_MT_density_pool(1:50:end),VIF_density_pool(1:50:end),'.');
%     title(['Movie-',num2str(iM),', scrabled MT and VIF densities Cross-correlation:', ...
%         num2str(corr(scrable_MT_density_pool,VIF_density_pool)),'(radius: ',num2str(radius),')']);
%     saveas(h13, [ML_ROOT_DIR,filesep,'scrable_MT_density_correlation_movie',num2str(iM),'(radius',num2str(radius),').tif']);
%     corr_array_S_MT(iM) = corr(scrable_MT_density_pool,VIF_density_pool);
%     
%     h14= figure(14);hold off;
%     plot(scrable_MT_density_pool(1:50:end),scrable_VIF_density_pool(1:50:end),'.');
%     title(['Movie-',num2str(iM),', scrabled MT and VIF densities Cross-correlation:', ...
%         num2str(corr(scrable_MT_density_pool,scrable_VIF_density_pool)),'(radius: ',num2str(radius),')']);
%     saveas(h14, [ML_ROOT_DIR,filesep,'scrable_MT_VIF_density_correlation_movie',num2str(iM),'(radius',num2str(radius),').tif']);
%     corr_array_SS(iM) = corr(scrable_MT_density_pool,scrable_VIF_density_pool);
%         
%     MT_density_all_pool = [MT_density_all_pool; MT_density_pool];
%     VIF_density_all_pool = [VIF_density_all_pool; VIF_density_pool];
%     
%     scrable_MT_density_all_pool = [scrable_MT_density_all_pool; scrable_MT_density_pool];
%     scrable_VIF_density_all_pool = [scrable_VIF_density_all_pool; scrable_VIF_density_pool];
%     
% end
% 
% h21= figure(21);
% plot(MT_density_all_pool(1:100:end),VIF_density_all_pool(1:100:end),'.');
% title(['MT and VIF densities Cross-correlation: ', ...
%     num2str(corr(MT_density_all_pool,VIF_density_all_pool)),'(radius: ',num2str(radius),')']);
% saveas(h21, [ML_ROOT_DIR,filesep,'scrable_VIF_density_correlation_movie',num2str(iM),'(radius',num2str(radius),').tif']);
% AA = axis;
% 
% h22= figure(22);hold off;
% plot(MT_density_all_pool(1:100:end),scrable_VIF_density_all_pool(1:100:end),'.');
% title(['Movie-',num2str(iM),', MT and scrabled VIF densities Cross-correlation:', ...
%     num2str(corr(MT_density_all_pool,scrable_VIF_density_all_pool)),'(radius: ',num2str(radius),')']);
% axis(AA);
% saveas(h22, [ML_ROOT_DIR,filesep,'scrable_VIF_density_correlation_movie',num2str(iM),'(radius',num2str(radius),').tif']);
% 
% h23= figure(23);hold off;
% plot(scrable_MT_density_all_pool(1:100:end),VIF_density_all_pool(1:100:end),'.');
% title(['Movie-',num2str(iM),', scrabled MT and VIF densities Cross-correlation:', ...
%     num2str(corr(scrable_MT_density_all_pool,VIF_density_all_pool)),'(radius: ',num2str(radius),')']);
% axis(AA);
% saveas(h23, [ML_ROOT_DIR,filesep,'scrable_MT_density_correlation_movie',num2str(iM),'(radius',num2str(radius),').tif']);
% 
% 
% h24= figure(24);hold off;
% plot(scrable_MT_density_all_pool(1:100:end),scrable_VIF_density_all_pool(1:100:end),'.');
% title(['Movie-',num2str(iM),', scrabled MT and VIF densities Cross-correlation:', ...
%     num2str(corr(scrable_MT_density_all_pool,scrable_VIF_density_all_pool)),'(radius: ',num2str(radius),')']);
% axis(AA);
% saveas(h24, [ML_ROOT_DIR,filesep,'scrable_MT_VIF_density_correlation_movie',num2str(iM),'(radius',num2str(radius),').tif']);
% 
% corr_ML = [corr(MT_density_all_pool,VIF_density_all_pool);
%     corr(MT_density_all_pool,scrable_VIF_density_all_pool);
%     corr(scrable_MT_density_all_pool,VIF_density_all_pool);
%     corr(scrable_MT_density_all_pool,scrable_VIF_density_all_pool);];
% 
% save([ML_ROOT_DIR,filesep,'ML_corr_radius_',num2str(radius),'.mat'],...
%     'corr_array','corr_array_S_VIF',...
%     'corr_array_S_MT','corr_array_SS',...
%     'corr_ML');

%%
