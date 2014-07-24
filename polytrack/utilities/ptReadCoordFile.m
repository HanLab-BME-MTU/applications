function coord = ptReadCoordFile (filename)
% ptReadCoordFile reads all coordinates from a file. The file should be in the format:
%         x1 y1\n
%         x2 y2\n
%         xn yn\n
%
% SYNOPSIS       coord = ptReadCoordFile (filename)
%
% INPUT          filename : the name of the file with the coordinates
%
% OUTPUT         coord: a matrix with size (n,2) where n is the number of coordinate pairs
%
% DEPENDENCIES   ptReadCoordFile uses { nothing }
%                                  
%                ptReadCoordFile is used by { }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Andre Kerstens        May 04          Initial Version

% Check that the file exists
if ~exist (filename, 'file')
   error ('The file with name %s does not exist!', filename);
   return;
end

% Open filename for reading
fid = fopen (filename, 'r');

% Read in the coordinates and transpose to get them in the format we want
coord = fscanf (fid,'%g %g',[2 inf]);
coord = coord';
