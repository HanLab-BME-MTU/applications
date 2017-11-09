function [deskewFileName]=unzipAndDeskewLatticeTimePoint(filename,outputDir)


system(sprintf('bunzip2 "%s"',filename));
[p,n,e]=fileparts(filename);
unzipFilename=[p filesep n];
vol=stackRead(unzipFilename);
[data_deskew] = shear3DinDim2(vol, 32.8, 0, 0.4, .1, 0, 0);
deskewFileName=[outputDir filesep n]; 
stackWrite(data_deskew,deskewFileName);