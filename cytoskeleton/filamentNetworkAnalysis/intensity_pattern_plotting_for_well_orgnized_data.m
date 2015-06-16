function [ML_pattern_median,ML_pattern_mean] = intensity_pattern_plotting_for_well_orgnized_data(ML, start_row, end_row, start_col, end_col)
% % after organized with each well and saved in ML, get intensity pattern

ML_ROOT_DIR = ML.outputDirectory_;

if(~exist('start_row','var'))
    start_row = 'A';
end

if(~exist('end_row','var'))
    end_row = 'P';
end

if(~exist('start_col','var'))
    start_col = 1;
end

if(~exist('end_col','var'))
    end_col = 24;
end


% the number of movies
movieNumber =  length(ML.movieDataFile_);

start_row_ind = start_row - 'A' + 1 ;
end_row_ind = end_row - 'A' + 1 ;

load(ML.movieDataFile_{1});
nChannel = numel(MD.channels_);

intensity_all_pool = cell(1,nChannel);

for iMD  = 1 : 8: movieNumber
    try
        % load this movie
        load(ML.movieDataFile_{iMD});
    catch
        disp('The MD file is missing');
        continue;
    end

    nFrame = MD.nFrames_;

    for iChannel = 1 : nChannel
        display(['Checking: iMD:', num2str(iMD), ', iChannel:', num2str(iChannel)]);

            for iFrame = 1 : nFrame
                display(['iMD:', num2str(iMD), ', iChannel:', num2str(iChannel),', iFrame:', num2str(iFrame)]);
                current_im  = MD.channels_(iChannel).loadImage(iFrame);
                small_img = current_im(1:20:end,1:20:end);
                intensity_all_pool{iChannel} = [intensity_all_pool{iChannel};small_img(:);];
            end% end of a frame

    end % end of a Channel

end % end of if previous gathering exists for MD

ML_mean_percentiles = cell(1, nChannel);
ML_pattern_median = cell(1, nChannel);
ML_pattern_mean = cell(1, nChannel);

% get the threshold for background
for iChannel = 1 : nChannel
    T_array(iChannel) = thresholdRosin(intensity_all_pool{iChannel});
    ML_mean_percentiles{iChannel} = nan(movieNumber,8);
    ML_pattern_median{iChannel} = nan(end_row_ind,end_col);
    ML_pattern_mean{iChannel} = nan(end_row_ind,end_col);
end
% save the theresholds
save([ML_ROOT_DIR,filesep,'intensity_stat_threshs.mat'],'T_array');


intensity_thisMovie_pool=cell(1,nChannel);


for iMD  = 1 : movieNumber
    try
        % load this movie
        load(ML.movieDataFile_{iMD});
    catch
        disp('The MD file is missing');
        continue;
    end
    
     MD_Dir = MD.outputDirectory_;
     
     if( MD_Dir(end)=='/' || MD_Dir(end)=='\' )
         MD_Dir = MD_Dir(1:end-1);
     end
     
    sind = find(MD_Dir=='\' | MD_Dir=='/');
    
    MD_Well_ID = MD_Dir(sind(end)+1:end);
    
    row_ind =  MD_Well_ID(6) - 'A' + 1;
    col_ind =  str2num(MD_Well_ID(8:9));
    
    if(row_ind>end_row_ind || row_ind< start_row_ind)
        continue;
    end
       
    if(col_ind>end_col || col_ind< start_col)
        continue;
    end
    
    nFrame = MD.nFrames_;
    
    for iChannel = 1 : nChannel
        intensity_thisMovie_pool{iChannel} = [];
        display(['Checking: iMD:', num2str(iMD), ', iChannel:', num2str(iChannel)]);
        
        for iFrame = 1 : nFrame
            display(['iMD:', num2str(iMD), ', iChannel:', num2str(iChannel),', iFrame:', num2str(iFrame)]);
            current_im  = MD.channels_(iChannel).loadImage(iFrame);
            intensity_thisMovie_pool{iChannel}  = [intensity_thisMovie_pool{iChannel};current_im(current_im>T_array(iChannel)); ];
        end % end of a frame
    ML_mean_percentiles{iChannel}(iMD,1:8) = mean_percentiles(intensity_thisMovie_pool{iChannel});    
    
    
    ML_pattern_median{iChannel}(row_ind, col_ind) = ML_mean_percentiles{iChannel}(iMD,5);
    ML_pattern_mean{iChannel}(row_ind, col_ind) = ML_mean_percentiles{iChannel}(iMD,1);
    
     end % end of a Channel    
    
end

save([ML_ROOT_DIR,filesep,'intensity_stat.mat'],'ML_pattern_median','ML_pattern_mean','T_array');


for iChannel = 1 : nChannel    
    
    h4=figure(4);
    plot_color_dot_from_matrix(ML_pattern_median{iChannel},25);
    
    set(gca,'XTick', (start_col:end_col)-0.5);
    
    for iW = 1 :end_col - start_col+1
        X_Tick_cell{iW} = num2str(iW +start_col-1);
    end
    set(gca,'XTickLabel',X_Tick_cell);
    
    set(gca,'YTick',(start_row_ind:end_row_ind)-0.5);
    for iM = 1 : end_row_ind - start_row_ind+1
        Y_Tick_cell{iM} = char('A'+iM + start_row_ind-1-1);
    end
    set(gca,'YTickLabel',Y_Tick_cell);
    
    plot([0 0],[end_row_ind round(end_row_ind*2)],'w','linewidth',10);
    plot([end_col round(end_col*2)],[0 0],'w','linewidth',10);
    
    saveas(h4,[ML_ROOT_DIR,filesep,'IntensityPattern_',num2str(iChannel),'-perwell.tif']);
    saveas(h4,[ML_ROOT_DIR,filesep,'IntensityPattern_',num2str(iChannel),'-perwell.fig']);
end


