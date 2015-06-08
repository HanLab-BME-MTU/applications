function amiraWriteTracks(filename,tracks,varargin)
% Write an Amira Mesh file with name [<filename>_%04d.am] representing tracks. 
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
ip.addParamValue('vertexProp',{}, @iscell);
ip.parse( varargin{:});
p=ip.Results;

[pathstr,name,ext] = fileparts(filename); 
basename=[pathstr filesep name];

s=ip.Results.scales;

parfor fIdx=1:tracks.numTimePoints
    
    %% Indx of tracks on the current frame
    tracksOn=([tracks.endFrame]>=fIdx)&(fIdx>[tracks.startFrame]);
    nbTracsOn=sum(tracksOn);
    
    %% tracks extremity
    tracksEnds=size(nbTracsOn*2,3);
    % relative Idx of the end of each tracks (e.g. for use in tracks(i).x)
    endRelIdx=fIdx-[tracks.startFrame]+1; 
    count=1;
    for tIdx=find(tracksOn)
        tr=tracks(tIdx);
        tracksEnds(count+1,1)=tr.x(endRelIdx(tIdx))*s(1);
        tracksEnds(count+1,2)=tr.y(endRelIdx(tIdx))*s(2);
        tracksEnds(count+1,3)=tr.z(endRelIdx(tIdx))*s(3);        
        tracksEnds(count,1)=tr.x(1)*s(1);
        tracksEnds(count,2)=tr.y(1)*s(2);
        tracksEnds(count,3)=tr.z(1)*s(3);                        
        count=count+2;
    end
    
    %% tracks extremity property (start or end)
    startEnd=ones(nbTracsOn*2,1);
    startEnd(1:2:end)=0;

    %% tracks edgeConnectivity
    tracksEdge=reshape(1:2*nbTracsOn,[2,nbTracsOn])'-1;    
    
    %% tracks point    
    % numPointsPerEdge=max(2,endRelIdx(tracksOn)'); %if one single time point Amira still needs two time points...
    numPointsPerEdge=endRelIdx(tracksOn)';
    numPoints=sum(numPointsPerEdge);
    tracksPoints=zeros(size(numPoints,3));
    count=1;
    for tIdx=find(tracksOn)
        tr=tracks(tIdx);
        endIdx=endRelIdx(tIdx);
        tracksPoints(count:(count+endIdx-1),1)=tr.x(1:endIdx)*s(1);
        tracksPoints(count:(count+endIdx-1),2)=tr.y(1:endIdx)*s(2);
        tracksPoints(count:(count+endIdx-1),3)=tr.z(1:endIdx)*s(3);        
        count=count+endIdx;
    end
    
    %% Track id (edge property)
    tracksId=find(tracksOn)'
    
    
    
    % Write 
    frameFilename=[basename '_t_' num2str(fIdx,'%04.0f'),'.am'];
    fid = fopen(frameFilename, 'w');
    fprintf(fid,['# Amira 2.0 ASCII\n\n']);
    fprintf(fid,['define VERTEX ' num2str(nbTracsOn*2) '\n']);
    fprintf(fid,['define EDGE ' num2str(nbTracsOn) ' \n']);
    fprintf(fid,['define POINT ' num2str(numPoints) '\n']);
    fprintf(fid,'\n');
    fprintf(fid,'Parameters { ContentType "HxSpatialGraph"}\n\n');
    fprintf(fid,'VERTEX { float[3] VertexCoordinates } @1\n');
    fprintf(fid,'EDGE { int[2] EdgeConnectivity } @2\n');
    fprintf(fid,'EDGE { int NumEdgePoints } @3\n');
    fprintf(fid,'POINT { float[3] EdgePointCoordinates } @4\n');
    fprintf(fid,'EDGE { int trackId } @5\n');
    fprintf(fid,'VERTEX {int startEnd } @6\n');
    fclose(fid);
    if(nbTracsOn)
        fid = fopen(frameFilename, 'a');        
        fprintf(fid,'\n@1\n');
        fclose(fid);
        dlmwrite(frameFilename, tracksEnds, '-append', 'delimiter',' ','precision', 16);
        
        fid = fopen(frameFilename, 'a');
        fprintf(fid,'\n@2\n');
        fclose(fid);
        dlmwrite(frameFilename, tracksEdge, '-append', 'delimiter',' ','precision', 16);
        
        fid = fopen(frameFilename, 'a');
        fprintf(fid,'\n@3\n');
        fclose(fid);
        dlmwrite(frameFilename, numPointsPerEdge, '-append', 'delimiter',' ','precision', 16);
        
        fid = fopen(frameFilename, 'a');
        fprintf(fid,'\n@4\n');
        fclose(fid);
        dlmwrite(frameFilename, tracksPoints, '-append', 'delimiter',' ','precision', 16);
        
        fid = fopen(frameFilename, 'a');
        fprintf(fid,'\n@5\n');
        fclose(fid);
        dlmwrite(frameFilename, tracksId, '-append', 'delimiter',' ','precision', 16);
        
        fid = fopen(frameFilename, 'a');
        fprintf(fid,'\n@6\n');
        fclose(fid);
        dlmwrite(frameFilename, startEnd, '-append', 'delimiter',' ','precision', 16);
        
        for propIdx=1:length(p.vertexProp)
            fid = fopen(frameFilename, 'a');
            fprintf(fid,['\n VERTEX { float ' p.vertexProp{propIdx}{1} ' } @' num2str(4+propIdx) '\n']);
            fprintf(fid,['@' num2str(4+propIdx) '\n']);
            fclose(fid);
            dlmwrite(frameFilename, p.vertexProp{propIdx}{2}{fIdx}, '-append', 'delimiter',' ','precision', 16)
        end
    end
end
    
    
