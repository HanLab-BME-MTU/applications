% This program reads the computed track file and generates the
% corresponding MATLAB data structure.
%
% Created by Ge Yang, LCCB, TSRI
%
% Version: 0.1
%
% Date: April 28, 2004
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Step 1: Read in the data file and generate all the data structure
function [track] = trackRead(filename)

fid = fopen(filename,'rt');

track_num = 0;
nonghost_track_num = 0;

max_track_len = 400; % Maximum number of frames is 400
track = struct('startID',0, 'len', 0, 'lifespan', 0, 'points', zeros(max_track_len, 2));

temp_array = [];
total_pts = 0;

while(~feof(fid))
    linetemp = fgetl(fid);
    linetemp_copy = linetemp; %save a copy

    % find out how many points are there
    indx1 = find(linetemp == '(');
    indx2 = find(linetemp == ')');
    point_num = length(indx1);
    
    % Clean out the markers
    linetemp(indx1) = ' ';
    linetemp(indx2) = ' ';
    linetemp(end) = ' ';
    
    temp_array = [temp_array; sscanf(linetemp,'%f')];
    
    if (linetemp_copy(end) == '.') % Now we have a complete track.
        %fprintf('A new track is now read in completely.\n');
        track_num = track_num + 1;
        track(track_num).startID = temp_array(1);
        track(track_num).len = temp_array(2);
       
        % A ghost track is a track with length one.
        if (temp_array(2) > 1)
            nonghost_track_num = nonghost_track_num + 1;
        end
        
        
        temp_array_len = length(temp_array);
        
        track(track_num).points = reshape(temp_array(3:temp_array_len) , 2, round(length(temp_array)/ 2 - 1))';
        % Remember that the first pair does NOT contain point coordinates.
        
%         [temp1, temp2] = size(track(track_num).points);
%         if (track(track_num).len ~= temp1)
%             error('Data integrity is violated.\n');
%         end
%         
        total_pts = total_pts + temp_array(2);
        if (floor(mod(total_pts, 5000)) == 0)
            fprintf('Another 5000 points have been read in.\n');
        %    fprintf('Number of points read in = %d', total_pts);
        %    pause;
        end
        temp_array = [];
    end
end

fprintf('Total number of pts is %d\n', total_pts);
fprintf('Total number of tracks is %d. Among them, the number of non-ghost tracks is %d\n', track_num, nonghost_track_num);

