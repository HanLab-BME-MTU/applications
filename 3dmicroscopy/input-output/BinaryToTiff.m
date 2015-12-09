%%
cd C:\Users\ALSM_Master\Desktop\Data\Reto\TestSpeed\0_75ms\151208\Cell9

fileID = fopen('1_CAM01_000001.bin','r+','l');

test=fread(fileID,'uint16');

fclose(fileID);

numRows=512;

numCols=128;

numZ=168;


test2=zeros(numRows,numCols,numZ);
test2=reshape(test,[numRows,numCols,numZ]);

figure;imshow(squeeze(sum(test2,2)),[])
figure;imshow(squeeze(max(permute(test2,[2,1,3]))),[])

bl=squeeze(max(permute(test2,[2,1,3])));