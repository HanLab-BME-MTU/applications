%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%small program to convert binary files to tiff
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
cd C:\Users\ALSM_Master\Desktop\Data\Reto\RPE_1ms\tractinGFP\151208\Cell10

fileID = fopen('1_CAM02_000002.bin','r+','l');

test=fread(fileID,'uint16');

fclose(fileID);

numRows=512;

numCols=64;

numZ=68;


test2=zeros(numRows,numCols,numZ);


test2=reshape(test,[numRows,numCols,numZ]);

%figure;imshow(squeeze(sum(test2,2)),[])
figure;imshow(squeeze(max(permute(test2,[2,1,3]))),[])

bl=squeeze(max(permute(test2,[2,1,3])));