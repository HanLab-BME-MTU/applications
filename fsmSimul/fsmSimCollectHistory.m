function stack = fsmSimCollectHistory(dirName,writeTiffs)
%FSMSIMCOLLECTHISTORY screens a directory after running fsmSimAssembly
%
% SYNOPSIS stack = fsmSimCollectHistory(dirName,writeTiffs)
% 
% INPUT dirName : directory name, if [] the screening is done in the current directory
%       writeTiffs : 0 if not, integer number <=16 to specify bitdepth
% 
% OUTPUT stack : data structure with the following fields 
%                (each field has the time as the last dimension)
%
%             .img : image stack 
%             .mesh: mesh with label distribution
%             .t   : timepoint
%             .N   : number of speckles
%             .sMap: map speckles

if isempty(dirName)
   oldDir = [];
else
   oldDir = cd(dirName);
end;

dirListing = dir;

% get file indices with the body 'fsmSimMaps'
iEntry = 1;
indxList = [];
for( i = 1:length(dirListing))
   if(~dirListing(i).isdir)
      if strncmp(dirListing(i).name,'fsmSimMaps',length('fsmSimMaps'))
         [path,body,no,ext]=getFilenameBody(dirListing(i).name);
         indxList(iEntry) = str2num(no);
         iEntry = iEntry + 1;
      end;
   end;
end;
indxList = sort(indxList);

for i = 1:length(indxList)
   load(sprintf('fsmSimMaps%.6d',indxList(i)));
   stack.t(i) = timePoint;
   stack.N(i) = nb;
   stack.img(:,:,i) = img;
   stack.lblMap(:,:,i) = lblMap;
   stack.sMap(:,:,i) = speckMap;
end;


if writeTiffs
   % determine bit depth
   if writeTiffs <=8
      imgStack = uint8(stack.img*2^writeTiffs);
   else
      if writeTiffs <=16
         imgStack = uint16(stack.img*2^writeTiffs);
      else
         error('writeTiffs exceeds bounds [0, 16]');
      end
   end;
   for i = 1:length(indxList)
      imwrite(imgStack(:,:,i),sprintf('fsmSimImg%.6d.tif',indxList(i)));
   end;
end;



      