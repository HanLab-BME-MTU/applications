function amiraWriteMovieInfo(filename, movieInfo,varargin)
% Write an Amira Mesh file with name [<filename>_%04d.am] representing vertex. 
% Options
%    - <scales>: [x y z] defines relative pixel size (must be synced to amira stack
%    opening)
%    - <prop>: {{'name',{Nx1,..}} ...} is the associate properties (Id,
%    probability ...)
%    
ip=inputParser();
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addParamValue('scales', [1 1 1], @isnumeric);
ip.addParamValue('prop',{}, @iscell);
ip.parse( varargin{:});
p=ip.Results;

[pathstr,name,ext] = fileparts(filename); 
basename=[pathstr name];
    
parfor fIdx=1:length(movieInfo)
    fMI=movieInfo(fIdx);
    numVertices=size(fMI.xCoord,1);
    %Write end points to Amira Surf object
    frameFilename=[filename '_t_' num2str(fIdx,'%04.0f'),'.am'];
    fid = fopen(frameFilename, 'w');
    fprintf(fid,['# Amira 2.0 ASCII\n\n']);
    fprintf(fid,['define VERTEX ',num2str(numVertices),'\n']);
    fprintf(fid,'define EDGE 1\n');
    fprintf(fid,'define POINT 2\n');
    fprintf(fid,'\n');
    fprintf(fid,'Parameters { ContentType "HxSpatialGraph"}\n\n');
    fprintf(fid,'EDGE { int[2] EdgeConnectivity } @1\n');
    fprintf(fid,'EDGE { int NumEdgePoints } @2\n');
    fprintf(fid,'POINT { float[3] EdgePointCoordinates } @3\n');
    fprintf(fid,'VERTEX { float[3] VertexCoordinates } @4\n');
    fprintf(fid,'\n');
    fprintf(fid,'@1\n');
    fprintf(fid,'0 0\n\n');
    fprintf(fid,'@2\n');
    fprintf(fid,'2\n\n');
    fprintf(fid,'@3\n');
    fprintf(fid,[num2str(fMI.xCoord(1,1)*p.scales(1)) ' ' num2str(fMI.yCoord(1,1)*p.scales(2)) ' ' num2str(fMI.zCoord(1,1)*p.scales(3)) '\n']);
    fprintf(fid,[num2str(fMI.xCoord(1,1)*p.scales(1)) ' ' num2str(fMI.yCoord(1,1)*p.scales(2)) ' ' num2str(fMI.zCoord(1,1)*p.scales(3)) '\n\n']);
    fprintf(fid,'@4\n');
    fclose(fid);
    dlmwrite(frameFilename, [fMI.xCoord(:,1)*p.scales(1) fMI.yCoord(:,1)*p.scales(2) fMI.zCoord(:,1)*p.scales(3)], '-append', 'delimiter',' ','precision', 16);
    for propIdx=1:length(p.prop)
        fid = fopen(frameFilename, 'a');
        fprintf(fid,['\n VERTEX { float ' p.prop{propIdx}{1} ' } @' num2str(4+propIdx) '\n']);
        fprintf(fid,['@' num2str(4+propIdx) '\n']);
        fclose(fid);
        dlmwrite(frameFilename, p.prop{propIdx}{2}{fIdx}, '-append', 'delimiter',' ','precision', 16)
    end
end 
    
    
