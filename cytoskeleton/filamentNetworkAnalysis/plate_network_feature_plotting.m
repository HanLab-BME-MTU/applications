% function plate_network_feature_plotting(Identifier_thisML,...
%     CFMP_feature_ordered_thisML,...
%     CFMP_feature_thisML,...
%     ChMP_feature_thisML,ChMP_feature_perwell_thisML, output_dir)
nFrame = size(CFMP_feature_ordered_thisML{1}{1},3);
iChannel = 1;
nMovie = numel(ChMP_feature_thisML);

ChMP_feature_thisML_ch1 = nan(18,8,nMovie);

for iMD = 1 : numel(ChMP_feature_thisML)
    ChMP_feature_thisML_ch1(:,:,iMD) = ChMP_feature_thisML{iMD}{1}(:,:);
end

whole_plate_stat = median(ChMP_feature_thisML_ch1,3);
whole_plate_stat(1:18,1) = mean(ChMP_feature_thisML_ch1(:,1,:),3);


for iMD = 1 : nMovie
    for iW = 1 : nFrame/4
        for iF = 1:18
            for iS = 1 :8
                CFMP_feature_ordered_thisCh_perwell(iMD,iW,iF,iS) = ...
                    median(CFMP_feature_ordered_thisML{iMD}{1}(iF,iS,(iW-1)*4+1:(iW)*4));
            end
        end
    end
end

for iMD = 1 : nMovie
    for iFrame = 1 : nFrame
        for iF = 1:18
            for iS = 1 :8
                CFMP_feature_ordered_thisCh_perframe(iMD,iFrame,iF,iS) = ...
                    median(CFMP_feature_ordered_thisML{iMD}{1}(iF,iS,iFrame));
            end
        end
    end
end

for iF = 1:18
    well_median = CFMP_feature_ordered_thisCh_perwell(:,:,iF,5);
    
    CFMP_feature_median_wells(iF) = prctile(well_median(:),50);
    
    CFMP_feature_5p_wells(iF) = prctile(well_median(:),5);
    
    CFMP_feature_95p_wells(iF) = prctile(well_median(:),95);
end

for iF = 1:18
    frame_median = CFMP_feature_ordered_thisCh_perframe(:,:,iF,5);
    
    CFMP_feature_median_frames(iF) = prctile(frame_median(:),50);
    
    CFMP_feature_5p_frames(iF) = prctile(frame_median(:),2);
    
    CFMP_feature_95p_frames(iF) = prctile(frame_median(:),98);
end
  

textcell = {'Straightness of filament (for each filament)',...
    'Length (for each filament)',...
    'Pixel number of segmentation (for each filament)',...
    'Filament density (for each valid pixel)',...
    'Scrabled filament density (for each valid pixel)',...
    'Filament orientation (for each pixel)',...
    'Filament orientation (for each pixel,centered)',...
    'Filamet intensity (integrated for each filament)',...
    'Filamet intensity (average for each filament)',...
    'Filamet intensity (integrated for each dilated filament)',...
    'Filamet intensity (average for each dilated filament)',...
    'Scale of detected filaments (for each pixel)',...
    'ST Responce (integrated for each filament)',...
    'ST Responce (average for each filament)',...
    'ST Responce (integrated for each dilated filament)',...
    'ST Responce (average for each dilated filament)',...
    'Filament curvature (average for each filament)',...
    'Filament curvature (for each pixel on filaments)'};

for iF =  [ 1 2 8:16] 
    % for the background, plot 25%~50%~75% patch
    % then the mean
%     h1 = figure(1); close;
%     h1 = figure(1); hold off;
%     p = patch([1  nMovie*nFrame/4 nMovie*nFrame/4 1 1],...
%         [whole_plate_stat(iF,4) whole_plate_stat(iF,4) whole_plate_stat(iF,6)...
%         whole_plate_stat(iF,6) whole_plate_stat(iF,4)],[0.8 0.9 1],...
%         'EdgeColor','none','linestyle','none');
%     
%     hold on;
%     %             line1 = plot([1 nMovie*nFrame/4],[whole_plate_stat(iF,1) whole_plate_stat(iF,1)],...
%     %                 'm','linewidth',2);
%     line2 = plot([1 nMovie*nFrame/4],[whole_plate_stat(iF,4) whole_plate_stat(iF,4)],...
%         '-.','color',[0.5 0.5 1]/1.02,'linewidth',2);
%     line3 = plot([1 nMovie*nFrame/4],[whole_plate_stat(iF,5) whole_plate_stat(iF,5)],...
%         '-','color',[0.5 0.5 1]/1.2,'linewidth',2);
%     line4 = plot([1 nMovie*nFrame/4],[whole_plate_stat(iF,6) whole_plate_stat(iF,6)],...
%         ':','color',[0.5 0.5 1]/1.02,'linewidth',2);
%     
%     %             legend([line1, line2,line3,line4],'Mean','25%','Median','50%');
%     legend([line2,line3,line4],'25%','Median','50%');
%     title(textcell{iF},'fontsize',15);
%     
    
%% First look at per site outliers

     h1 = figure(1); hold off; plot(0,0);
    h1 = figure(1); hold off;
    p = patch([1  nMovie*nFrame nMovie*nFrame 1 1],...
        [CFMP_feature_5p_frames(iF) CFMP_feature_5p_frames(iF) CFMP_feature_95p_frames(iF)...
        CFMP_feature_95p_frames(iF) CFMP_feature_5p_frames(iF)],[0.8 0.9 1],...
        'EdgeColor','none','linestyle','none');
    
    hold on;
    %             line1 = plot([1 nMovie*nFrame],[whole_plate_stat(iF,1) whole_plate_stat(iF,1)],...
    %                 'm','linewidth',2);
    line2 = plot([1 nMovie*nFrame],[CFMP_feature_5p_frames(iF) CFMP_feature_5p_frames(iF)],...
        '-.','color',[0.5 0.5 1]/1.02,'linewidth',2);
    line3 = plot([1 nMovie*nFrame],[CFMP_feature_median_frames(iF) CFMP_feature_median_frames(iF)],...
        '-','color',[0.5 0.5 1]/1.2,'linewidth',2);
    line4 = plot([1 nMovie*nFrame],[CFMP_feature_95p_frames(iF) CFMP_feature_95p_frames(iF)],...
        ':','color',[0.5 0.5 1]/1.02,'linewidth',2);
     %             legend([line1, line2,line3,line4],'Mean','25%','Median','50%');
    legend([line2,line3,line4],'2%','Median','98%');
    title(textcell{iF},'fontsize',15);
    
    
    % for samples
    % plot median for each one
    
    Y = CFMP_feature_ordered_thisCh_perframe(:,:,iF,5);
    X = 1 :numel(Y);
    Y = Y';
    Y = Y(:);
    X_in = X(Y>= CFMP_feature_5p_frames(iF) & Y<= CFMP_feature_95p_frames(iF));
    Y_in = Y(Y>= CFMP_feature_5p_frames(iF) & Y<= CFMP_feature_95p_frames(iF));
    X_out = X(Y< CFMP_feature_5p_frames(iF) | Y> CFMP_feature_95p_frames(iF));
    Y_out = Y(Y< CFMP_feature_5p_frames(iF) | Y> CFMP_feature_95p_frames(iF));
        
    
    plot(X_in,Y_in,'*');
    plot(X_out, Y_out,'r*');
    axis([0 max(X)+1 ...
        min(CFMP_feature_median_frames(iF)+1.1*(min(Y)-CFMP_feature_median_frames(iF)),CFMP_feature_median_frames(iF)+2.5*(CFMP_feature_5p_frames(iF)-CFMP_feature_median_frames(iF)))...
        max(CFMP_feature_median_frames(iF)+1.1*(max(Y)-CFMP_feature_median_frames(iF)),CFMP_feature_median_frames(iF)+2.5*(CFMP_feature_95p_frames(iF)-CFMP_feature_median_frames(iF)))]);
    
    Identifier_together = cat(1,Identifier_thisML{:});
    
    for iX = 1:numel(X_out)
        iXX = X_out(iX);
        [iframe,im] = ind2sub([nFrame, nMovie,], iXX);
        ID  = Identifier_together{im, iframe};
        ID = ID(end-4:end-0);
        Row_charatcter = char('A'+ (ID(1)-'0')*10+ID(2)-'0'+1);    
        
        if(Y_out(iX)>CFMP_feature_median_wells(iF))
        text(iXX-0.8, Y_out(iX)+(max(Y)-min(Y))*0.03, [Row_charatcter,'-',ID(3:4),'-',ID(5)],'color',[1 0 0]);
        else
        text(iXX-0.8, Y_out(iX)-(max(Y)-min(Y))*0.03, [Row_charatcter,'-',ID(3:4),'-',ID(5)],'color',[1 0 0]);
        end
    end
saveas(h1,[output_dir,filesep,'Feature_',num2str(iF),'-persite.jpg']);
    saveas(h1,[output_dir,filesep,'Feature_',num2str(iF),'-persite.tiff']);
    saveas(h1,[output_dir,filesep,'Feature_',num2str(iF),'-persite.fig']);
    

%% Second the per well outliers (without deleting the frame outliers)

     h2 = figure(2); hold off; plot(0,0);
    h2 = figure(2); hold off;
    p = patch([1  nMovie*nFrame/4 nMovie*nFrame/4 1 1],...
        [CFMP_feature_5p_wells(iF) CFMP_feature_5p_wells(iF) CFMP_feature_95p_wells(iF)...
        CFMP_feature_95p_wells(iF) CFMP_feature_5p_wells(iF)],[0.8 0.9 1],...
        'EdgeColor','none','linestyle','none');
    
    hold on;
    %             line1 = plot([1 nMovie*nFrame/4],[whole_plate_stat(iF,1) whole_plate_stat(iF,1)],...
    %                 'm','linewidth',2);
    line2 = plot([1 nMovie*nFrame/4],[CFMP_feature_5p_wells(iF) CFMP_feature_5p_wells(iF)],...
        '-.','color',[0.5 0.5 1]/1.02,'linewidth',2);
    line3 = plot([1 nMovie*nFrame/4],[CFMP_feature_median_wells(iF) CFMP_feature_median_wells(iF)],...
        '-','color',[0.5 0.5 1]/1.2,'linewidth',2);
    line4 = plot([1 nMovie*nFrame/4],[CFMP_feature_95p_wells(iF) CFMP_feature_95p_wells(iF)],...
        ':','color',[0.5 0.5 1]/1.02,'linewidth',2);
     %             legend([line1, line2,line3,line4],'Mean','25%','Median','50%');
    legend([line2,line3,line4],'5%','Median','95%');
    title(textcell{iF},'fontsize',15);
    
    
    % for samples
    % plot median for each one
    
    Y = CFMP_feature_ordered_thisCh_perwell(:,:,iF,5);
    Y=Y';
    X = 1 :numel(Y);
    Y = Y(:);
    X_in = X(Y>= CFMP_feature_5p_wells(iF) & Y<= CFMP_feature_95p_wells(iF));
    Y_in = Y(Y>= CFMP_feature_5p_wells(iF) & Y<= CFMP_feature_95p_wells(iF));
    X_out = X(Y< CFMP_feature_5p_wells(iF) | Y> CFMP_feature_95p_wells(iF));
    Y_out = Y(Y< CFMP_feature_5p_wells(iF) | Y> CFMP_feature_95p_wells(iF));
    
    
    
    plot(X_in,Y_in,'*');
    plot(X_out, Y_out,'r*');
    axis([0 max(X)+1 ...
        min(CFMP_feature_median_wells(iF)+1.1*(min(Y)-CFMP_feature_median_wells(iF)),CFMP_feature_median_wells(iF)+2.5*(CFMP_feature_5p_wells(iF)-CFMP_feature_median_wells(iF)))...
        max(CFMP_feature_median_wells(iF)+1.1*(max(Y)-CFMP_feature_median_wells(iF)),CFMP_feature_median_wells(iF)+2.5*(CFMP_feature_95p_wells(iF)-CFMP_feature_median_wells(iF)))]);
    
    Identifier_together = cat(1,Identifier_thisML{:});
    
    for iX = 1:numel(X_out)
        iXX = X_out(iX);
        [iw,im] = ind2sub([nFrame/4,nMovie], iXX);
        ID  = Identifier_together{im, iw*4};
        ID = ID(end-4:end-1);
        Row_charatcter = char('A'+ (ID(1)-'0')*10+ID(2)-'0'+1);    
        
        if(Y_out(iX)>CFMP_feature_median_wells(iF))
        text(iXX-0.8, Y_out(iX)+(max(Y)-min(Y))*0.03, [Row_charatcter,'-',ID(3:4)],'color',[1 0 0]);
        else
        text(iXX-0.8, Y_out(iX)-(max(Y)-min(Y))*0.03, [Row_charatcter,'-',ID(3:4)],'color',[1 0 0]);
        end
    end
    saveas(h2,[output_dir,filesep,'Feature_',num2str(iF),'-perwell.jpg']);
    saveas(h2,[output_dir,filesep,'Feature_',num2str(iF),'-perwell.tiff']);
    saveas(h2,[output_dir,filesep,'Feature_',num2str(iF),'-perwell.fig']);
    
    h4=figure(4);
    plot_color_dot_from_matrix(squeeze(...
        CFMP_feature_ordered_thisCh_perframe(:,:,iF,5)),25);
    
    set(gca,'XTick',2:4:nFrame);
    for iW = 1 :nFrame/4
    X_Tick_cell{iW} = num2str(iW);
    end
    set(gca,'XTickLabel',X_Tick_cell);

     set(gca,'YTick',0.5 :nMovie-0.5);
    for iM = 1 :nMovie
    Y_Tick_cell{iM} = char('A'+iM-1);
    end
    set(gca,'YTickLabel',Y_Tick_cell);
    
    for iFrame = 0:4:nFrame
     plot([iFrame iFrame],[0 nMovie]);
    end
    plot([0 0],[nMovie round(nMovie*2)],'w','linewidth',10);
    plot([nFrame round(nFrame*2)],[0 0],'w','linewidth',10);
    
    saveas(h4,[output_dir,filesep,'FeatureDots_',num2str(iF),'-persite.jpg']);
    saveas(h4,[output_dir,filesep,'FeatureDots_',num2str(iF),'-persite.tiff']);
    saveas(h4,[output_dir,filesep,'FeatureDots_',num2str(iF),'-persite.fig']);
    
        h5=figure(5);
    plot_color_dot_from_matrix(squeeze(...
        CFMP_feature_ordered_thisCh_perwell(:,:,iF,5)),25);
    
    set(gca,'XTick',0.5:nFrame/4-0.5);
    for iW = 1 :nFrame/4
    X_Tick_cell{iW} = num2str(iW);
    end
    set(gca,'XTickLabel',X_Tick_cell);

     set(gca,'YTick',0.5 :nMovie-0.5);
    for iM = 1 :nMovie
    Y_Tick_cell{iM} = char('A'+iM-1);
    end
    set(gca,'YTickLabel',Y_Tick_cell);
    plot([0 0],[nMovie round(nMovie*2)],'w','linewidth',10);
    plot([nFrame/4 round(nFrame)],[0 0],'w','linewidth',10);
   
    saveas(h5,[output_dir,filesep,'FeatureDots_',num2str(iF),'-perwell.jpg']);
    saveas(h5,[output_dir,filesep,'FeatureDots_',num2str(iF),'-perwell.tiff']);
    saveas(h5,[output_dir,filesep,'FeatureDots_',num2str(iF),'-perwell.fig']);
  

end