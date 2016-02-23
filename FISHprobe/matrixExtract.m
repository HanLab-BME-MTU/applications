function [Ch1,Ch2,Ch3] = matrixExtract(movieName)


Path = strcat('', movieName);
MD = MovieData.load(Path);
xSize = MD.imSize_(1);
ySize = MD.imSize_(2);
nDepth = MD.zSize_;
frameN = 1;

for Cha=1:3
    zStack3D = zeros(xSize, ySize, nDepth);
    for i = 1:nDepth
        zStack3D(:,:,i) = MD.getChannel(Cha).loadImage(frameN,i);
    end
    zStack3D=100*zStack3D/max(zStack3D(:));
    
    switch Cha
        case 1
            Ch1=zStack3D;
        case 2
            Ch2=zStack3D;
        case 3
            Ch3=zStack3D;
    end           
end