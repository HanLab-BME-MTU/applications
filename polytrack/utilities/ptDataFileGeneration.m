% program to generate data file to be processed by the LAP tracker.

startID = 25;
endID = 50;
frame_num = endID - startID + 1;

frames = struct('points', [], 'length', 0);
for i = startID: endID
    filename = strcat('coordinates0', num2str(i), '.txt');
    frames(i-startID+1).points = load(filename);
    frames(i-startID+1).length = length(frames(i-startID+1).points);
end

datafile_name = 'poly_001_junk.txt';

fprintf('Generating speckle image file %s\n', filename);
fid = fopen(datafile_name, 'wt');
total_pt_num = 0;

h = waitbar(0, 'Generating point data files, Please wait...');
for i = 1 : frame_num
        total_pt_num = total_pt_num + frames(i).length;
        
        for j = 1 : frames(i).length
            temp = round(frames(i).points(j, 1:2));
            fprintf(fid,' (');
            fprintf(fid, '%10.6f %10.6f', temp(2), temp(1));
            
            if ((temp(2) < 0) | (temp(1)<0))
                fprintf('Negative coordinates. %f  %f\n', temp(2), temp(1));
                pause;
                pause;
                pause;
            end
            fprintf(fid,')');
            if (mod(j,5) == 0)
                fprintf(fid,' \n');
            end
        end
        fprintf(fid,' .\n'); 
        waitbar(i / frame_num, h);
end
fclose(fid);
close(h);
fprintf('Frame point file %s has been generated.\n', datafile_name);
fprintf('Total number of points is %d\n', total_pt_num);
