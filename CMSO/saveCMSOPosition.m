function saveCMSOPosition(csvFilename,movieInfo,timeInterval,channel)

fid = fopen(csvFilename, 'wt');
fprintf(fid,'#COM An average detection of 3D intracellular object detection.\n');
fprintf(fid,'#ONT OME http://www.openmicroscopy.org/Schemas/Documentation/Generated/OME-2015-01/ome.html\n');
fprintf(fid,'#ONT CMSO http://cmso.org/...\n');
fprintf(fid,['CMSO:ElapsedMS,\t' 'OME:X,\t' 'OME:Y,\t' 'OME:Z,\t' 'OME:C,\t' 'CMSO:Intensity' '\n']);
 
for fIdx=1:length(movieInfo)
    if(isfield(movieInfo(fIdx),'zCoord'))
        Z=movieInfo(fIdx).zCoord(:,1);
    else
        Z=zeros(length(movieInfo(fIdx).xCoord(:,1)),1);
    end
    for row = 1:length(movieInfo(fIdx).xCoord(:,1))
      fprintf(fid, '%0.2f,\t%0.2f,\t%0.2f,\t%0.2f,\t%d,\t%0.2f\n',timeInterval*fIdx, movieInfo(fIdx).xCoord(row,1), movieInfo(fIdx).yCoord(row,1),Z(row),channel,movieInfo(fIdx).int(row,1));
    end
end

fclose(fid);

