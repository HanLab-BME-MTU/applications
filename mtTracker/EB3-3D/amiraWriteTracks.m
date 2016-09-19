function amiraWriteTracks(filename,tracks,varargin)
% Write an Amira Mesh file with name [<filename>_%04d.am] representing tracks. 
% Options
%    - <scales>: [x y z t] defines relative pixel size (must be synced to
%    amira stack opening)
%    - <movieData>: alternative to scale input (WARNING: the scale option
%    has priority)
%    - <vertexProp>: {{'name',{NVertex x 1, ...}},...} is the vertex-associated 
%    properties  each frame must be described
%    in a cell. 
%    - <edgeProp>: {{'name',NTrack x 1}, ...} is the
%    edge-associated properties
%    
ip=inputParser();
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addParamValue('scales', [1 1 1 1], @isnumeric);
ip.addParamValue('MD',[],@(MD) isa(MD,'MovieData'));
ip.addParamValue('vertexProp',{}, @iscell);
ip.addParamValue('fillGaps',true, @islogical);
ip.addParamValue('edgeProp',{}, @iscell);
ip.parse( varargin{:});
p=ip.Results;

[pathstr,name,ext] = fileparts(filename); 
basename=[pathstr filesep name];


s=ip.Results.scales;
if(length(s)==3) s=[s 1]; end; %BW comp check. 

if ((all(s==[1 1 1 1]))&&(~isempty(p.MD)))
    s=[p.MD.pixelSize_ p.MD.pixelSize_ p.MD.pixelSizeZ_ p.MD.timeInterval_];
end

if(~exist(fileparts(filename))) mkdir(fileparts(filename)); end;

% GAP filling using the last known position (gap are still mark by tracksFeatIndxCG
% trackFeat
if(p.fillGaps)
    se=[zeros(1,tracks.numTimePoints) 1 ones(1,tracks.numTimePoints)];
    for tIdx=1:length(tracks)
        gi=tracks(tIdx).gapMask;
        if(any(gi))
            copyIdx=1:tracks(tIdx).lifetime;
            copyIdx(gi)=0;
            copyIdx=imdilate(copyIdx,se);
            tracks(tIdx).x=tracks(tIdx).x(copyIdx);
            tracks(tIdx).y=tracks(tIdx).y(copyIdx);
            tracks(tIdx).z=tracks(tIdx).z(copyIdx);
        end
    end
end

% Frame 0 is the cumulative track distribution. Ugly ? I know, stfu.
for fIdx=0:tracks.numTimePoints
    
    %% Indx of tracks on the current frame
    if(fIdx>0)
        tracksOn=([tracks.endFrame]>=fIdx)&(fIdx>[tracks.startFrame]);
    else % 0 frames is the cumulative display
        tracksOn=true(1,length(tracks));
    end
    nbTracsOn=sum(tracksOn);
    
    %% tracks extremity
    tracksEnds=size(nbTracsOn*2,3);
    % relative Idx of the end of each tracks (e.g. for use in tracks(i).x)
    endRelIdx=fIdx-[tracks.startFrame]+1; 
    if(fIdx==0)
        endRelIdx=[tracks.lifetime];
    end;
    count=1;
    for tIdx=find(tracksOn)
        tr=tracks(tIdx);
        tracksEnds(count+1,1)=(tr.x(endRelIdx(tIdx))-1)*s(1);
        tracksEnds(count+1,2)=(tr.y(endRelIdx(tIdx))-1)*s(2);
        tracksEnds(count+1,3)=(tr.z(endRelIdx(tIdx))-1)*s(3);
        tracksEnds(count,1)=(tr.x(1)-1)*s(1);
        tracksEnds(count,2)=(tr.y(1)-1)*s(2);
        tracksEnds(count,3)=(tr.z(1)-1)*s(3);                        
        count=count+2;
    end

    
    %% tracks extremity property (start or end)
    startEnd=ones(nbTracsOn*2,1);
    startEnd(1:2:end)=0;

    %% tracks edgeConnectivity
    tracksEdge=reshape(1:2*nbTracsOn,[2,nbTracsOn])'-1;    
    
    %% tracks point    
    % numPointsPerEdge=max(2,endRelIdx(tracksOn)'); %if one single time point Amira still needs two time points...
    numPointsPerEdge=endRelIdx((tracksOn))';
    numPoints=sum(numPointsPerEdge);
    tracksPoints=zeros([numPoints,3]);
    pointType=zeros(numPoints,1);
    count=1;
    for tIdx=find(tracksOn)
        tr=tracks(tIdx);
        endIdx=endRelIdx(tIdx);
        gapM=tr.gapMask;
        tracksPoints(count:(count+endIdx-1),1)=(tr.x(1:endIdx)-1)*s(1);
        tracksPoints(count:(count+endIdx-1),2)=(tr.y(1:endIdx)-1)*s(2);
        tracksPoints(count:(count+endIdx-1),3)=(tr.z(1:endIdx)-1)*s(3);      
        pointType(count:(count+endIdx-1))=gapM(1:endIdx)';
        count=count+endIdx;
    end
    
    %% Track id (edge property)
    tracksId=find(tracksOn)';
    
    %% Track lifetime (edge property)
    tracksLft=[tracks(tracksOn).lifetime]';
    
    %% Track Median Speed (edge property)
    tracksMedSpeed= arrayfun(@(t) nanmedian(sum((   [s(1)*t.x(1:end-1);s(2)*t.y(1:end-1);s(3)*t.z(1:end-1)]- ... 
                                                    [s(1)*t.x(2:end);  s(2)*t.y(2:end); s(3)*t.z(2:end)  ]).^2).^0.5/s(4)) ,tracks(tracksOn));
    %% Track Max Speed (edge property)
    tracksMaxSpeed= arrayfun(@(t)    nanmax(sum((   [s(1)*t.x(1:end-1);s(2)*t.y(1:end-1);s(3)*t.z(1:end-1)]- ... 
                                                    [s(1)*t.x(2:end);  s(2)*t.y(2:end); s(3)*t.z(2:end)  ]).^2).^0.5/s(4)) ,tracks(tracksOn),'unif',0);
                                                
    %% Track diffCoeff (edge property)
    tracksDiffCoeff=arrayfun(@(t) nanmean(sum([s(1)*t.x(1)-s(1)*t.x(2:end); s(2)*t.y(1)-s(2)*t.y(2:end); s(3)*t.z(1)-s(3)*t.z(2:end)].^2))/(6*t.lifetime*s(4)) ,tracks(tracksOn));
    
    % Write 
    paramCount=0;
    frameFilename=[basename '_t_' num2str(fIdx,'%04.0f'),'.am'];
    if(fIdx==0)
        frameFilename=[basename '_t_' num2str(fIdx,'%04.0f'), '_cumulative_display.am'];
    end
    fid = fopen(frameFilename, 'w');
    fprintf(fid,['# AmiraMesh 3D ASCII 2.0\n\n']);
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
    fprintf(fid,'EDGE { int lifetime } @7\n');    
    fprintf(fid,'EDGE { float medianSpeed} @8\n');  
    fprintf(fid,'EDGE { float maxSpeed} @9\n');  
    fprintf(fid,'EDGE { float diffCoeff} @10\n');
    fprintf(fid,'POINT { int pointType } @11\n');
    fclose(fid);
    if(nbTracsOn)
        paramCount=paramCount+1;
        fid = fopen(frameFilename, 'a');        
        fprintf(fid,'\n@1\n');
        fclose(fid);
        dlmwrite(frameFilename, tracksEnds, '-append', 'delimiter',' ','precision', 16);
        
        paramCount=paramCount+1;
        fid = fopen(frameFilename, 'a');
        fprintf(fid,'\n@2\n');
        fclose(fid);
        dlmwrite(frameFilename, tracksEdge, '-append', 'delimiter',' ','precision', 16);
        
        paramCount=paramCount+1;
        fid = fopen(frameFilename, 'a');
        fprintf(fid,'\n@3\n');
        fclose(fid);
        dlmwrite(frameFilename, numPointsPerEdge, '-append', 'delimiter',' ','precision', 16);
        
        paramCount=paramCount+1;
        fid = fopen(frameFilename, 'a');
        fprintf(fid,'\n@4\n');
        fclose(fid);
        dlmwrite(frameFilename, tracksPoints, '-append', 'delimiter',' ','precision', 16);
        
        paramCount=paramCount+1;
        fid = fopen(frameFilename, 'a');
        fprintf(fid,'\n@5\n');
        fclose(fid);
        dlmwrite(frameFilename, tracksId, '-append', 'delimiter',' ','precision', 16);
        
        paramCount=paramCount+1;
        fid = fopen(frameFilename, 'a');
        fprintf(fid,'\n@6\n');
        fclose(fid);
        dlmwrite(frameFilename, startEnd, '-append', 'delimiter',' ','precision', 16);

        paramCount=paramCount+1;
        fid = fopen(frameFilename, 'a');
        fprintf(fid,'\n@7\n');
        fclose(fid);
        dlmwrite(frameFilename, tracksLft, '-append', 'delimiter',' ','precision', 16);
        
        paramCount=paramCount+1;
        fid = fopen(frameFilename, 'a');
        fprintf(fid,'\n@8\n');
        fclose(fid);
        dlmwrite(frameFilename, tracksMedSpeed, '-append', 'delimiter',' ','precision', 16);

        paramCount=paramCount+1;
        fid = fopen(frameFilename, 'a');
        fprintf(fid,'\n@9\n');
        fclose(fid);
        dlmwrite(frameFilename, tracksMaxSpeed, '-append', 'delimiter',' ','precision', 16);
        
        paramCount=paramCount+1;
        fid = fopen(frameFilename, 'a');
        fprintf(fid,'\n@10\n');
        fclose(fid);
        dlmwrite(frameFilename, tracksDiffCoeff, '-append', 'delimiter',' ','precision', 16);
        
        paramCount=paramCount+1;
        fid = fopen(frameFilename, 'a');
        fprintf(fid,'\n@11\n');
        fclose(fid);
        dlmwrite(frameFilename, pointType, '-append', 'delimiter',' ','precision', 16);
                
        for propIdx=1:length(p.vertexProp)
            fid = fopen(frameFilename, 'a');
            paramCount=paramCount+1;
            fprintf(fid,['\nVERTEX { float ' p.vertexProp{propIdx}{1} ' } @' num2str(paramCount+propIdx-1) '\n']);
            fprintf(fid,['@' num2str(paramCount+propIdx-1) '\n']);
            fclose(fid);
            dlmwrite(frameFilename, p.vertexProp{propIdx}{2}{fIdx}, '-append', 'delimiter',' ','precision', 16)           
        end
        
        for propIdx=1:length(p.edgeProp)
            fid = fopen(frameFilename, 'a');
            paramCount=paramCount+1;
            fprintf(fid,['\nEDGE { float ' p.edgeProp{propIdx}{1} ' } @' num2str(paramCount+propIdx-1) '\n']);
            fprintf(fid,['@' num2str(paramCount+propIdx-1) '\n']);
            fclose(fid);
            dlmwrite(frameFilename, p.edgeProp{propIdx}{2}(tracksOn), '-append', 'delimiter',' ','precision', 16)
        end
    end
end
    
    
