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

if(~exist(fileparts(filename))) mkdir(fileparts(filename)); end;


[pathstr,name,ext] = fileparts(filename); 
basename=[pathstr filesep name];

parfor fIdx=1:length(movieInfo)
    fMI=movieInfo(fIdx);
    numVertices=size(fMI.xCoord,1);
    %Write end points to Amira Surf object
    headerString={};
    headerString=[headerString;{'# Amira 2.0 ASCII '}];
    headerString=[headerString;{['define VERTEX ',num2str(numVertices),' ']}];
    headerString=[headerString;'define EDGE 1 '];
    headerString=[headerString;'define POINT 2 '];
    headerString=[headerString;' '];
    headerString=[headerString;'Parameters { ContentType "HxSpatialGraph"} '];
    headerString=[headerString;'EDGE { int[2] EdgeConnectivity } @1 '];
    headerString=[headerString;'EDGE { int NumEdgePoints } @2'];
    headerString=[headerString;'POINT { float[3] EdgePointCoordinates } @3'];
    headerString=[headerString;'VERTEX { float[3] VertexCoordinates } @4'];
    headerString=[headerString;' '];
    headerString=[headerString;'@1'];
    headerString=[headerString;'0 0'];
    headerString=[headerString;' '];
    headerString=[headerString;'@2'];
    headerString=[headerString;'2'];
    headerString=[headerString;' '];
    headerString=[headerString;'@3'];
    
    detectionString={};
    if(~isempty(fMI.xCoord))
        %dlmwrite(frameFilename, repmat([(fMI.xCoord(1,1)-1)*p.scales(1) (fMI.yCoord(1,1)-1)*p.scales(2) (fMI.zCoord(1,1)-1)*p.scales(3)],2,1), '-append', 'delimiter',' ','precision', 16);
        detectionString={strjoin(cellstr(num2str(repmat([(fMI.xCoord(1,1)-1)*p.scales(1) (fMI.yCoord(1,1)-1)*p.scales(2) (fMI.zCoord(1,1)-1)*p.scales(3)],2,1)))','\n')};
    end
%     fid = fopen(frameFilename, 'a');
%     headerSpring=[headerSpring;'@4\n');
%     fclose(fid);
    detectionString=[detectionString; ' '; '@4'];
    if(~isempty(fMI.xCoord))
       %dlmwrite(frameFilename, [(fMI.xCoord(:,1)-1)*p.scales(1) (fMI.yCoord(:,1)-1)*p.scales(2) (fMI.zCoord(:,1)-1)*p.scales(3)], '-append', 'delimiter',' ','precision', 16);
       detectionString=[detectionString; strjoin(cellstr(num2str([(fMI.xCoord(:,1)-1)*p.scales(1) (fMI.yCoord(:,1)-1)*p.scales(2) (fMI.zCoord(:,1)-1)*p.scales(3)]))','\n')];
    end
    propString={};
    for propIdx=1:length(p.prop)
%         fid = fopen(frameFilename, 'a');
%         fprintf(fid,['\n VERTEX { float ' p.prop{propIdx}{1} ' } @' num2str(4+propIdx) '\n']);
%         fprintf(fid,['@' num2str(4+propIdx) '\n']);
%         fclose(fid);
%         dlmwrite(frameFilename, p.prop{propIdx}{2}{fIdx}, '-append', 'delimiter',' ','precision', 16)
        
        propString=[propString; '\n VERTEX { float ' p.prop{propIdx}{1} ' } @' num2str(4+propIdx) '\n'];
        propString=[propString; ['@' num2str(4+propIdx) '\n']];
        propString=[propString; strjoin(cellstr(num2str([(fMI.xCoord(:,1)-1)*p.scales(1) (fMI.yCoord(:,1)-1)*p.scales(2) (fMI.zCoord(:,1)-1)*p.scales(3)]))','\n')];
    end
    frameFilename=[basename '_t_' num2str(fIdx,'%04.0f'),'.am'];
    fid = fopen(frameFilename, 'w');
    fprintf(fid,'%s\n', strjoin(headerString','\n'));
    fprintf(fid,'%s\n', strjoin(detectionString','\n'));
    fprintf(fid,'%s\n', strjoin(propString','\n'));
    fclose(fid);
end 
    
% cumulFrameFilename=[basename '_cumulative.am'];
% cfid = fopen(cumulFrameFilename, 'w');
% fprintf(cfid,['# Amira 2.0 ASCII\n\n']);
% fprintf(cfid,['define VERTEX ',num2str(numVertices),'\n']);
% fprintf(cfid,'define EDGE 1\n');
% fprintf(cfid,'define POINT 2\n');
% fprintf(cfid,'\n');
% fprintf(cfid,'Parameters { ContentType "HxSpatialGraph"}\n\n');
% fprintf(cfid,'EDGE { int[2] EdgeConnectivity } @1\n');
% fprintf(cfid,'EDGE { int NumEdgePoints } @2\n');
% fprintf(cfid,'POINT { float[3] EdgePointCoordinates } @3\n');
% fprintf(cfid,'VERTEX { float[3] VertexCoordinates } @4\n');
% fprintf(cfid,'\n');
% fprintf(cfid,'@1\n');
% fprintf(cfid,'0 0\n\n');
% fprintf(cfid,'@2\n');
% fprintf(cfid,'2\n\n');
% fprintf(cfid,'\n@3\n');


    
