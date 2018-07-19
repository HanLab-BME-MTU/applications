function [deskewFileName]=deskewLatticeTimePoint(filename,outputDir)
[p,n,e]=fileparts(filename);
vol=stackRead(filename);
[data_deskew] = shear3DinDim2(vol, 32.8, 0, 0.4, .1, 0, 0);
deskewFileName=[outputDir filesep n e]; 
stackWrite(data_deskew,deskewFileName);