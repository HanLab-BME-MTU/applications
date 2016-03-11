%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%small program to convert binary files to 3D Stack
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function im=BinaryTo3D(fileName,Xdim,Ydim,Zdim)


fileID = fopen(fileName);

test=fread(fileID,'uint16');
size(test)
fclose(fileID);


im=zeros(Xdim,Ydim,Zdim);


im=reshape(test,[Xdim,Ydim,Zdim]);

%figure;imshow(squeeze(sum(test2,2)),[])
%figure;imshow(squeeze(max(permute(im,[2,1,3]))),[])

