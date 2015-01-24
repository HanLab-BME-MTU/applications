function amiraWriteMovieInfo(filename, movieInfo)
% Write an Amira Mesh file with name [<filename>_%04d.am] representing vertex. 

[pathstr,name,ext] = fileparts(filename); 
basename=[pathstr name];
    
parfor fIdx=1:length(movieInfo)
    fMI=movieInfo(fIdx);
    numVertices=size(fMI.xCoord,1);
    %Write end points to Amira Surf object
    frameFilename=[filename '_t_' num2str(fIdx,'%04.0f'),'.am'];
    fid = fopen(frameFilename, 'w');
    fprintf(fid,['# HyperSurface 0.1 ASCII\n\n Vertices ',num2str(numVertices),'\n']);
    fprintf(fid,['\n']);
    fclose(fid);
    dlmwrite(frameFilename, [fMI.xCoord(:,1) fMI.yCoord(:,1) fMI.zCoord(:,1)], '-append', 'delimiter',' ');
end 
    
    
