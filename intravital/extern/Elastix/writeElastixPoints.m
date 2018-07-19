function writeElastixPoints(points, filename, format)
% WRITEELASTIXTRANSFORM Write Elastix transforms parameters in text files.
% Author: Paul Balança
%
% WRITEELASTIXTRANSFORM(elastixTransforms, filename)
%
%     Input:
%        points              Matrix of points to write (size: Nx2 or Nx3).
%        filename            Filename where to write points file.
%        format              Format of the points: 'index' (in moving volume) or 'point' (in mm).
%

% Write points file
filePoints = fopen(filename, 'wt');
if filePoints == -1
    error('Could not open file "%s".\n', filename);
end

% Index : convert to C
if strcmp(format, 'index')
    points = points - 1;
end

% First two lines
fprintf(filePoints, '%s\n', format);
fprintf(filePoints, '%d\n', size(points, 1));

% Write lines (points coordinates)
for i = 1:size(points, 1)
    if size(points, 2) == 2
        fprintf(filePoints, '%.12f %.12f\n', points(i,1), points(i,2));
    else
        fprintf(filePoints, '%.12f %.12f %.12f\n', points(i,1), points(i,2), points(i,3));        
    end
end
fclose(filePoints);

end
