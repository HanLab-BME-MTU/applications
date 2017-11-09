% function  protrusion_hole_correction(MD)
% This function is to prevent protrusion spikes to form enclosed holes
% only limit one cell in the field


display_msg_flag = 0; % display warning or not
package_process_ind_script;

H2=fspecial('disk', 2);
H2 = H2>0;

H3=fspecial('disk', 2);
H3 = H3>0;

for iChannel = 1 : 1 %numel(MD.channels_);
    
    output_dir =  MD.processes_{indexCellRefineProcess}.outFilePaths_{iChannel};
    output_dir_clean = [output_dir,  '_clean'];
    output_dir_restored = [output_dir,  '_restored'];
    mkdir(output_dir_clean);
    mkdir(output_dir_restored);
    
    Channel_FilesNames = MD.channels_(iChannel).getImageFileNames(1:MD.nFrames_);
    
    for iFrame = 1 : MD.nFrames_
        
        Cell_Mask = ((MD.processes_{indexCellSegProcess}.loadChannelOutput(iChannel,iFrame))>0);
        mask = keep_largest_area(Cell_Mask);
        
        open_mask  = imopen(mask,H2);
        open_mask = keep_largest_area( open_mask);
        
%         if(sum(sum(imfill(open_mask,'holes') - open_mask))>0)
%           open_mask  = imopen(mask,H3);
%           open_mask = keep_largest_area( open_mask);
%         end
        
        mask_diff = (( mask) - ( open_mask));
        
        [diff_label, num] = bwlabel(mask_diff);
        
        restored_mask = open_mask;
        
        for n = 1 : num
            added_mask = or(restored_mask,diff_label==n);
            if(sum(sum(imfill(added_mask,'holes') - or(imfill(restored_mask,'holes'), added_mask)))==0)
                restored_mask = added_mask;
            end
        end
        
        imwrite(open_mask,[output_dir_clean,filesep,'refined_mask_',Channel_FilesNames{iFrame}]);
        imwrite(restored_mask,[output_dir_restored,filesep,'refined_mask_',Channel_FilesNames{iFrame}]);
        
    end
end
