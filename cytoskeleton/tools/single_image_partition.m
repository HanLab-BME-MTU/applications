function outImg_Cell = single_image_partition(inImg, row_partition_number, col_partition_number)
% Function of partitioning single image into multiple smaller ones, with
% user input
%           
% Input:      inImg:           the input image to be partitioned
%             row_partition_number:     how many parts in the row
%                                       direction, that is the y-axis;
%                                       default 3
%             col_partition_number:     how many parts in the column
%                                       direction, that is the x-axis;
%                                       default 3

% Output:     outImg_Cell:     a cell structure each cell contains the
%                              partitioned images

% Created 06 2014 by Liya Ding, Matlab R2012b

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
if(isempty(inImg))
    % if is empty, return empty
    outImg_Cell=[];
    return;
end

if(iscell(inImg))
     % if input is a cell, return empty
    outImg_Cell=[];
    return;
end

%% main function code
outImg_Cell = cell(row_partition_number, col_partition_number);


img_size_y = size(inImg,1);
img_size_x = size(inImg,2);

part_size_y = round(img_size_y/row_partition_number);

part_y_start = 1:part_size_y:img_size_y;

% correct for row number cannot be exactly divided
if(numel(part_y_start)>row_partition_number)
    % there could one more than requested, so get rid of the last one
    part_y_start = part_y_start(1:row_partition_number);
else
    if(part_y_start(row_partition_number)+part_size_y-1 > img_size_y)
    % in this case the last piece is smaller than others, so make it the
    % same with some overlap
    part_y_start(row_partition_number) = img_size_y - (part_size_y-1);
    end
end

%same for rows

part_size_x = round(img_size_x/col_partition_number);

part_x_start = 1:part_size_x:img_size_x;

% correct for col number cannot be exactlx divided
if(numel(part_x_start)>col_partition_number)
    % there could one more than requested, so get rid of the last one
    part_x_start = part_x_start(1:col_partition_number);
else
    if(part_x_start(col_partition_number)+part_size_x-1 > img_size_x)
    % in this case the last piece is smaller than others, so make it the
    % same with some overlap
    part_x_start(col_partition_number) = img_size_x - (part_size_x-1);
    end
end


% get the output cell 
for iR = 1 : row_partition_number
    for iC = 1 : col_partition_number
        outImg_Cell{iR, iC} = inImg(part_y_start(iR):part_y_start(iR)+part_size_y-1, ...
            part_x_start(iC):part_x_start(iC)+part_size_x-1);
    end
end

