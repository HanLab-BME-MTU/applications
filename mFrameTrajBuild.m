function MFT = mFrameTrajBuild(M,fpt)
%MFRAMETRAJBUILD  Create trajectories of speckles over multiframes.
%
% SYNOPSIS : 
%    MFT = mFrameSpeckTrack(M,fpt)
%
% INPUT :
%    M : The vector field of the displacements of speckles created by Aaron's 
%       Speckle Tracking Software. It is an m-by-4-by-n multi-dimentional array
%       where 'm' is roughly (due to the death of old speckles and the birth of 
%       new speckles) the number of speckles and 'n+1' is the total number of 
%       frames. The 4 columns have the form [y0 x0 y x] where (y0,x0) is the
%       base and (y,x) is the tip of the displacement vector.
%    fpt : The number of frames per trajectory. It is the number of frames over
%       which the trajectory of each speckle is tracked.
%
% OUTPUT :
%    MFT : A cell array of 'N' elements where 'N' equals the rounded integer
%       of the total number of frames divided by 'fpt'. The kth element of the
%       cell array is an m-by-2*fpt matrix where 'm' is the number of speckles
%       of frame 'fpt*(k-1)+1'. The '2*fpt' columns have the form 
%       [y1 x1 y2 x2 ... yn xn] where (yi,xi) is the tracked position of the 
%       speckle in frame 'fpt*(k-1)+i'.

%Get the total number of frames in 'M'.
numFrames = size(M,3)+1;

if fpt > numFrames | fpt < 2
   error(['The number of frames per trajectory, ''fpt'' has to be ' ...
      'great than 1 and less than the total number of frames stored in ' ...
      '''M'' that is size(M,3)+1.']);
end

N = floor(numFrames/fpt);

MFT = cell(1,N);

%The kth element of 'MFT' stores the trajectories of speckles from frame
% 'fpt*(k-1)+1' to frame 'fpt*k'.
corLen  = 2;
frameID = 1;
for k = 1:N
   %We need to use 'find' to remove some zeros due to the death of old speckles 
   % and the birth of new speckles.
   MFT{k}(:,1:4) = M(find(M(:,1,frameID)~=0 & M(:,3,frameID)~=0),:,frameID);
   jj = 3;
   for i = 3:fpt
      frameID = frameID+1;
      MFT{k}(:,jj:jj+3) = vectorFieldInterp( ...
         M(find(M(:,1,frameID)~=0 & M(:,3,frameID)~=0),:,frameID), ...
         MFT{k}(:,jj:jj+1),corLen,[]);
      jj = jj+2;
   end
   frameID = frameID+1;
end
