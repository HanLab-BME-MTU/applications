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

numVerticesCell=cell(1,length(movieInfo));
cumulDetectionValues=cell(1,length(movieInfo));
cumulPropValues=cell(length(p.prop)+1,length(movieInfo));
for fIdx=1:length(movieInfo)
    fMI=movieInfo(fIdx);
    numVertices=size(fMI.xCoord,1);
    numVerticesCell{fIdx}=numVertices;
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
    
    detectionValues={};
    detectionString={};
    if(~isempty(fMI.xCoord))
        detectionValues={strjoin(cellstr(num2str(repmat([(fMI.xCoord(1,1)-1)*p.scales(1) (fMI.yCoord(1,1)-1)*p.scales(2) (fMI.zCoord(1,1)-1)*p.scales(3)],2,1)))','\n')};
    end
    detectionString=[detectionValues; ' '; '@4'];

    detectionValues={};
    if(~isempty(fMI.xCoord))
       detectionValues=strjoin(cellstr(num2str([(fMI.xCoord(:,1)-1)*p.scales(1) (fMI.yCoord(:,1)-1)*p.scales(2) (fMI.zCoord(:,1)-1)*p.scales(3)]))','\n');
       detectionString=[detectionString;detectionValues];
    end
    cumulDetectionValues{fIdx}=detectionValues;

    propString={};
    propValues={};
    
    propIdx=1;
    propString=[propString; 'VERTEX {float vertexId} @5'];
    propString=[propString; ['@5']];
    propValues=strjoin(cellstr(num2str((1:numVertices)'))','\n');
    propString=[propString;propValues];
    cumulPropValues{1,fIdx}=propValues;
    
    for optPropIdx=1:length(p.prop)
        propIdx=propIdx+1;
        propString=[propString; 'VERTEX { float ' p.prop{optPropIdx}{1} ' } @' num2str(5+optPropIdx)'];
        propString=[propString; ['@' num2str(5+optPropIdx) ]];
        propValues=strjoin(cellstr(num2str([p.prop{optPropIdx}{2}{fIdx}]))','\n');
        propString=[propString;propValues];
        cumulPropValues{propIdx,fIdx}=propValues;
        
    end
    
    frameFilename=[basename '_t_' num2str(fIdx,'%04.0f'),'.am'];
    fid = fopen(frameFilename, 'w');
    fprintf(fid,'%s\n', strjoin(headerString','\n'));
    fprintf(fid,'%s\n', strjoin(detectionString','\n') );
    fprintf(fid,'%s\n', strjoin(propString','\n') );
    fclose(fid);
end 

headerString={};
headerString=[headerString;{'# Amira 2.0 ASCII '}];
headerString=[headerString;{['define VERTEX ',num2str(sum([numVerticesCell{:}])),' ']}];
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

detectionValues={};
detectionString={};
if(~isempty(fMI.xCoord))
    detectionValues={strjoin(cellstr(num2str(repmat([(movieInfo(1).xCoord(1,1)-1)*p.scales(1) (movieInfo(1).yCoord(1,1)-1)*p.scales(2) (movieInfo(1).zCoord(1,1)-1)*p.scales(3)],2,1)))','\n')};
end
detectionString=[detectionValues; ' '; '@4'];

detectionValues={};
if(~isempty(fMI.xCoord))
    detectionValues=strjoin(cumulDetectionValues,'\n');
    detectionString=[detectionString;detectionValues];
end

propString={};
propValues={};
for propIdx=1:length(p.prop)
    
    propString=[propString; 'VERTEX { float ' p.prop{propIdx}{1} ' } @' num2str(4+propIdx)'];
    propString=[propString; ['@' num2str(4+propIdx) ]];
    %propValues=strjoin(cumulPropValues{propIdx,:},'\n');
    propValues=[cumulPropValues{propIdx,:}];
    propString=[propString;propValues];
end

frameFilename=[basename '_cumulative.am'];
fid = fopen(frameFilename, 'w');
fprintf(fid,'%s\n', strjoin(headerString','\n'));
fprintf(fid,'%s\n', strjoin(detectionString','\n') );
fprintf(fid,'%s\n', strjoin(propString','\n') );
fclose(fid);

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


    
