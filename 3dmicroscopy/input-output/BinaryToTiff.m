%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%small program to convert binary files to tiff
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function im=BinaryToTiff(fileName,Xdim,Ydim,Zdim)


fileID = fopen(fileName);

test=fread(fileID,'uint16');

fclose(fileID);


im=zeros(Xdim,Ydim,Zdim);


im=reshape(test,[Xdim,Ydim,Zdim]);

%figure;imshow(squeeze(sum(test2,2)),[])
figure;imshow(squeeze(max(permute(im,[2,1,3]))),[])

