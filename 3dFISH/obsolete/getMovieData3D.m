function [MD,dataProperties] = getMovieData3D(moviePath, Cha)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Path = strcat('', movieName);
MD = MovieData.load(Path);
xSize = MD.imSize_(1);
ySize = MD.imSize_(2);
nDepth = MD.zSize_;
frameN = MD.nFrames_;

end

