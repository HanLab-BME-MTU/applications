function [] = TiffWrite(matrix,filename);

% June 1, 2003, by DK
% Saves 3D matrix as multi-image TIFF
% will overwrite existing tiff file without asking
% usage:  TiffWrite(matrix,filename);

if exist(filename) == 2
    
    delete(filename)
end


% for i=1:size(matrix,3)
%     imwrite(matrix(:,:,i),filename,'Compression','none','WriteMode','append')
% end

imwrite(matrix,filename, 'TIF','Compression','none','WriteMode','append','ColorSpace','rgb')
