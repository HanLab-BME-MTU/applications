function plot_differnce_map_wrapper(difference_map,outdir,customized_outdir,iFrame,radius,show_save_everything_flag,save_tif_flag)
% wrapper function for plotting out the network comparison details

if(~exist('save_tif_flag','var'))
    save_tif_flag = 0;
end


%%
% nan as white
difference_map.similarity_scoremap_proximity(difference_map.similarity_scoremap_proximity<0.1)=nan;
if(show_save_everything_flag==1)
    h3=figure(3);hold off;
    imagesc_white_nan(difference_map.similarity_scoremap_proximity,0,1);axis image;axis off;colorbar;
    CCC = colormap;
    
    colorbar;
    
    if CCC(1,1)> CCC(end,1)
        flip_colormap; display('Checkpoint 3');        
    end
    axis ij;
    title('Proximity Scoremap');
    %       outdir
    
    %%
    if(save_tif_flag>0)
        saveas(h3,[outdir,filesep,'proximity_score_frame_',num2str(iFrame),'tif']);
    end
    
    saveas(h3,[outdir,filesep,'proximity_score_frame_',num2str(iFrame),'.fig']);
    
    %%
    if(save_tif_flag>0)
        saveas(h3,[customized_outdir,filesep,'proximity_score_frame_',num2str(iFrame),'.tif']);
    end
    saveas(h3,[customized_outdir,filesep,'proximity_score_frame_',num2str(iFrame),'.fig']);
    
    %%
    difference_map.similarity_scoremap_alignment(difference_map.similarity_scoremap_alignment<0.1)=nan;
    h4=figure(4);hold off;
    imagesc_white_nan(difference_map.similarity_scoremap_alignment,0,1);axis image;axis off;colorbar;
    CCC = colormap;
    if CCC(1,1)> CCC(end,1)
        flip_colormap;
    end
    axis ij;
    title('Alignment Scoremap');
    if(save_tif_flag>0)           saveas(h4,[outdir,filesep,'alignment_score_frame_',num2str(iFrame),'.tif']); end
    saveas(h4,[outdir,filesep,'alignment_score_frame_',num2str(iFrame),'.fig']);
    if(save_tif_flag>0)     saveas(h4,[customized_outdir,filesep,'alignment_score_frame_',num2str(iFrame),'.tif']); end
    saveas(h4,[customized_outdir,filesep,'alignment_score_frame_',num2str(iFrame),'.fig']);
    %   display('4');
    difference_map.similarity_scoremap_combined((difference_map.similarity_scoremap_combined)<0.1)=nan;
   
    
    %%
    h5=figure(5);hold off;
    imagesc_white_nan(difference_map.similarity_scoremap_combined,0,1);axis image;axis off;colorbar;
    CCC = colormap;
    if CCC(1,1)> CCC(end,1)
        flip_colormap;
    end
    axis ij;
    title('Combined Scoremap');
    if(save_tif_flag>0)     saveas(h5,[outdir,filesep,'combined_score_frame_',num2str(iFrame),'.tif']); end
    saveas(h5,[outdir,filesep,'combined_score_frame_',num2str(iFrame),'.fig']);
    if(save_tif_flag>0)     saveas(h5,[customized_outdir,filesep,'combined_score_frame_',num2str(iFrame),'.tif']); end
    saveas(h5,[customized_outdir,filesep,'combined_score_frame_',num2str(iFrame),'.fig']);
    
end

%

if(show_save_everything_flag==1)
    h3=figure(3);
    imagesc_white_nan(difference_map.distance_map_1_2,0,radius*1);axis image;axis off;colorbar;
    flip_colormap;
    title('Distance Measure 1->2');
    axis ij;
    if(save_tif_flag>0)     saveas(h3,[outdir,filesep,'Dis12_frame_',num2str(iFrame),'.tif']); end
    saveas(h3,[outdir,filesep,'Dis12_frame_',num2str(iFrame),'.fig']);
    if(save_tif_flag>0)     saveas(h3,[customized_outdir,filesep,'Dis12_frame_',num2str(iFrame),'.tif']); end
    saveas(h3,[customized_outdir,filesep,'Dis12_frame_',num2str(iFrame),'.fig']);
    
    
    h3=figure(3);
    imagesc_white_nan(difference_map.distance_map_2_1,0,radius*1);axis image;axis off;colorbar;
    flip_colormap;
    title('Distance Measure 2->1');
    axis ij;
    if(save_tif_flag>0)     saveas(h3,[outdir,filesep,'Dis21_frame_',num2str(iFrame),'.tif']); end
    saveas(h3,[outdir,filesep,'Dis21_frame_',num2str(iFrame),'.fig']);
    
    if(save_tif_flag>0)     saveas(h3,[customized_outdir,filesep,'Dis21_frame_',num2str(iFrame),'.tif']); end
    saveas(h3,[customized_outdir,filesep,'Dis21_frame_',num2str(iFrame),'.fig']);
    
   
    h3=figure(3);
    imagesc_white_nan(difference_map.angle_map_2_1,0,pi/(2));axis image;axis off;colorbar;
    flip_colormap;
    title('orientation Measure 1->2');
    axis ij;
    if(save_tif_flag>0)     saveas(h3,[outdir,filesep,'Ang12_frame_',num2str(iFrame),'.tif']); end
    saveas(h3,[outdir,filesep,'Ang12_frame_',num2str(iFrame),'.fig']);
    if(save_tif_flag>0)     saveas(h3,[customized_outdir,filesep,'Ang12_frame_',num2str(iFrame),'.tif']); end
    saveas(h3,[customized_outdir,filesep,'Ang12_frame_',num2str(iFrame),'.fig']);
    
    h3=figure(3);
    imagesc_white_nan(difference_map.angle_map_1_2,0,pi/(2));axis image;axis off;colorbar;
    flip_colormap;
    title('orientation Measure 2->1');
    axis ij;
    if(save_tif_flag>0)     saveas(h3,[outdir,filesep,'Ang21_frame_',num2str(iFrame),'.tif']); end
    saveas(h3,[outdir,filesep,'Ang21_frame_',num2str(iFrame),'.fig']);
    if(save_tif_flag>0)     saveas(h3,[customized_outdir,filesep,'Ang21_frame_',num2str(iFrame),'.tif']); end
    saveas(h3,[customized_outdir,filesep,'Ang21_frame_',num2str(iFrame),'.fig']);
end


% show_angle12 = difference_map.angle_map_1_2;
% show_angle12(isnan(show_angle12)) = 30;
% show_angle21 = difference_map.angle_map_2_1;
% show_angle21(isnan(show_angle21)) = 30;
% show_dis12 = difference_map.distance_map_1_2;
% show_dis12(isnan(show_dis12)) = radius*1;
% show_dis21 = difference_map.distance_map_2_1;
% show_dis21(isnan(show_dis21)) = radius*1;

% if(show_save_everything_flag==1)
%     h3=figure(3);
%     imagesc(show_angle12+show_angle21+show_dis12+show_dis21);axis image;axis off;colorbar;
%     flip_colormap;
%     title('Sum of all Measures');
%  axis ij;
%     saveas(h3,[outdir,filesep,'AngDisSum_frame_',num2str(iFrame),'.tif']);
%     saveas(h3,[outdir,filesep,'AngDisSum_frame_',num2str(iFrame),'.fig']);
%    saveas(h3,[customized_outdir,filesep,'AngDisSum_frame_',num2str(iFrame),'.tif']);
%     saveas(h3,[customized_outdir,filesep,'AngDisSum_frame_',num2str(iFrame),'.fig']);
% end


% display the local supported distance and angle different matrix
if(show_save_everything_flag==1)
    h3=figure(3);
    imagesc_white_nan(difference_map.score_maps_distance_1_2,0,radius);axis image;axis off;colorbar;
    flip_colormap;
    title('Distance Measure 1->2 with Local Support');
    axis ij;
    if(save_tif_flag>0)     saveas(h3,[outdir,filesep,'Smooth_Dis12_frame_',num2str(iFrame),'.tif']); end
    saveas(h3,[outdir,filesep,'Smooth_Dis12_frame_',num2str(iFrame),'.fig']);
    if(save_tif_flag>0)     saveas(h3,[customized_outdir,filesep,'Smooth_Dis12_frame_',num2str(iFrame),'.tif']); end
    saveas(h3,[customized_outdir,filesep,'Smooth_Dis12_frame_',num2str(iFrame),'.fig']);
    
    h3=figure(3);
    imagesc_white_nan(difference_map.score_maps_distance_2_1,0,radius);axis image;axis off;colorbar;
    title('Distance Measure 2->1 with Local Support');
    axis ij;
    flip_colormap;
    if(save_tif_flag>0)     saveas(h3,[outdir,filesep,'Smooth_Dis21_frame_',num2str(iFrame),'.tif']); end
    saveas(h3,[outdir,filesep,'Smooth_Dis21_frame_',num2str(iFrame),'.fig']);
    if(save_tif_flag>0)      saveas(h3,[customized_outdir,filesep,'Smooth_Dis21_frame_',num2str(iFrame),'.tif']); end
    saveas(h3,[customized_outdir,filesep,'Smooth_Dis21_frame_',num2str(iFrame),'.fig']);
    
    h3=figure(3);
    imagesc_white_nan(difference_map.score_maps_angle_1_2,0,pi/2);axis image;axis off;colorbar;
    flip_colormap;
    title('Orientation Measure 1->2 with Local Support');
    axis ij;
    if(save_tif_flag>0)     saveas(h3,[outdir,filesep,'Smooth_Ang12_frame_',num2str(iFrame),'.tif']); end
    saveas(h3,[outdir,filesep,'Smooth_Ang12_frame_',num2str(iFrame),'.fig']);
    if(save_tif_flag>0)     saveas(h3,[customized_outdir,filesep,'Smooth_Ang12_frame_',num2str(iFrame),'.tif']); end
    saveas(h3,[customized_outdir,filesep,'Smooth_Ang12_frame_',num2str(iFrame),'.fig']);
    
    h3=figure(3);
    imagesc_white_nan(difference_map.score_maps_angle_2_1,0,pi/2);axis image;axis off;colorbar;
    flip_colormap;
    title('Orientation Measure 2->1 with Local Support');
    axis ij;
    if(save_tif_flag>0)     saveas(h3,[outdir,filesep,'Smooth_Ang21_frame_',num2str(iFrame),'.tif']);  end
    saveas(h3,[outdir,filesep,'Smooth_Ang21_frame_',num2str(iFrame),'.fig']);
    if(save_tif_flag>0)     saveas(h3,[customized_outdir,filesep,'Smooth_Ang21_frame_',num2str(iFrame),'.tif']); end
    saveas(h3,[customized_outdir,filesep,'Smooth_Ang21_frame_',num2str(iFrame),'.fig']);
    
    
    h4=figure(4); imagesc_white_nan((difference_map.score_maps_distance_2_1+difference_map.score_maps_distance_1_2)/2,0,radius);axis image;axis off;
    flip_colormap;
    title('Distance Measure 1->2 + 2->1 with Local Support');
    axis ij;colorbar;
    if(save_tif_flag>0)     saveas(h4,[outdir,filesep,'Final_D_frame_',num2str(iFrame),'.tif']); end
    saveas(h4,[outdir,filesep,'Final_D_frame_',num2str(iFrame),'.fig']);
    if(save_tif_flag>0)     saveas(h4,[customized_outdir,filesep,'Final_D_frame_',num2str(iFrame),'.tif']); end
    saveas(h4,[customized_outdir,filesep,'Final_D_frame_',num2str(iFrame),'.fig']);
    
    
    h5=figure(5); imagesc_white_nan(abs(difference_map.score_maps_angle_2_1/2)+abs(difference_map.score_maps_angle_1_2/2),0,1);axis image;axis off;
    flip_colormap;
    title('Orientation Measure 1->2 + 2->1 with Local Support');
    axis ij;colorbar;
    if(save_tif_flag>0)     saveas(h5,[outdir,filesep,'Final_A_frame_',num2str(iFrame),'.tif']);  end
    saveas(h5,[outdir,filesep,'Final_A_frame_',num2str(iFrame),'.fig']);
    if(save_tif_flag>0)     saveas(h5,[customized_outdir,filesep,'Final_A_frame_',num2str(iFrame),'.tif']); end
    saveas(h5,[customized_outdir,filesep,'Final_A_frame_',num2str(iFrame),'.fig']);
end



if(show_save_everything_flag==1)
    % similarity_scoremap(similarity_scoremap<0.2)=0.2;
    h6=figure(6); imagesc_white_nan(difference_map.similarity_scoremap_combined,0,1);
    axis image;axis off; axis ij;colorbar;
    title(['Similarity Score for frame ',num2str(iFrame)]);
    if(save_tif_flag>0)     saveas(h6,[outdir,filesep,'Network_sm_score_frame_',num2str(iFrame),'.tif']); end
    saveas(h6,[outdir,filesep,'Network_sm_score_frame_',num2str(iFrame),'.fig']);
    if(save_tif_flag>0)     saveas(h6,[customized_outdir,filesep,'Network_sm_score_frame_',num2str(iFrame),'.tif']); end
    saveas(h6,[customized_outdir,filesep,'Network_sm_score_frame_',num2str(iFrame),'.fig']);
    
    % similarity_scoremap(similarity_scoremap<0.2)=0.2;
    h7=figure(7); imagesc_white_nan(difference_map.similarity_scoremap_1to2,0,1);
    axis image;axis off; axis ij;colorbar;
    title(['Similarity Score 1to2 for frame ',num2str(iFrame)]);
    if(save_tif_flag>0)     saveas(h7,[outdir,filesep,'Network_1to2_sm_score_frame_',num2str(iFrame),'.tif']); end
    saveas(h7,[outdir,filesep,'Network_1to2_sm_score_frame_',num2str(iFrame),'.fig']);
    if(save_tif_flag>0)     saveas(h7,[customized_outdir,filesep,'Network_1to2_sm_score_frame_',num2str(iFrame),'.tif']); end
    saveas(7,[customized_outdir,filesep,'Network_1to2_sm_score_frame_',num2str(iFrame),'.fig']);
    
    % similarity_scoremap(similarity_scoremap<0.2)=0.2;
    h8=figure(8); imagesc_white_nan(difference_map.similarity_scoremap_2to1,0,1);
    axis image;axis off; axis ij;colorbar;
    title(['Similarity Score 2to1 for frame ',num2str(iFrame)]);
    if(save_tif_flag>0)     saveas(h8,[outdir,filesep,'Network_2to1_sm_score_frame_',num2str(iFrame),'.tif']); end
    saveas(h8,[outdir,filesep,'Network_2to1_sm_score_frame_',num2str(iFrame),'.fig']);
    if(save_tif_flag>0)     saveas(h8,[customized_outdir,filesep,'Network_2to1_sm_score_frame_',num2str(iFrame),'.tif']); end
    saveas(h8,[customized_outdir,filesep,'Network_2to1_sm_score_frame_',num2str(iFrame),'.fig']);
    
end
