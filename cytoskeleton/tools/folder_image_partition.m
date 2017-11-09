function folder_image_partition(image_folder, row_partition_number, col_partition_number, output_dir)
% Function of partitioning single image into multiple smaller ones, with
% user input
%           
% Input:      image_folder:           the input folder containing images to be partitioned
%             row_partition_number:     how many parts in the row
%                                       direction, that is the y-axis;
%                                       default 3
%             col_partition_number:     how many parts in the column
%                                       direction, that is the x-axis;
%                                       default 3
%             output_dir:               where you want the images to be
%                                       saved.

% Output:     None, all images are save to disk

% Created 07 2014 by Liya Ding, Matlab R2012b

%%
% if no input of row partition number, set it as default 3
if(nargin<2)
   row_partition_number=3;
end

% if no input of col partition number, set it as default 3
if(nargin<3)
   col_partition_number=3;
end

% correct less than 1 case and also float number case

row_partition_number = max(1, round(row_partition_number));
col_partition_number = max(1, round(col_partition_number));

% check input image
if(isempty(image_folder))
    % if is empty, return empty
    return;
end

if(iscell(image_folder))
     % if input is a cell, return empty
   return;
end

% if there is no defined output dir, have a default one
if(nargin<4)
   output_dir=image_folder;
   if(output_dir(end)==filesep)
       output_dir = output_dir(1:end-1);
   end
   output_dir = [output_dir, '_partitioned_images'];  
   
end

if(~exist(output_dir,'dir'))
    mkdir(output_dir);
end

% find all the images in this folder
input_dir_content1 = dir(fullfile([image_folder,filesep,'*.tif']));
input_dir_content2 = dir(fullfile([image_folder,filesep,'*.tiff']));
input_dir_content3 = dir(fullfile([image_folder,filesep,'*.jpg']));
input_dir_content4 = dir(fullfile([image_folder,filesep,'*.bmp']));

input_dir_content = [input_dir_content1;input_dir_content2;input_dir_content3;...
    input_dir_content4];

% the number of all images
nFrame = numel(input_dir_content);

for iFrame = 1 : nFrame
    inImg = imread([image_folder,filesep,input_dir_content(iFrame).name]);
    %     info = imfinfo([image_folder,filesep,input_dir_content(iFrame).name]);
    
    % do the image partition for this frame
    outImg_Cell = single_image_partition(inImg, row_partition_number, col_partition_number);
    
    for iR = 1 : row_partition_number
        for iC = 1 : col_partition_number
            % save each image partitioned for this frame
            imwrite(outImg_Cell{iR, iC},[output_dir,filesep, ...
                'image_frame_', num2str(iFrame),'_R_',num2str(iR), '_C_',num2str(iC), ...
                '.tif']);
        end
    end
    
    % % if there is anything related to bits later
    %     if(info.BitDepth == 16)
    %
    %     else
    %         if(info.BitDepth == 8)
    %         end
    %
    %     end
end


