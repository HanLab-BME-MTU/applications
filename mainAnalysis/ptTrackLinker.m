function [MPM,M]=ptTrackLinker(M)
% ptTrackLinker creates the magic position matrix MPM from M
%
% SYNOPSIS      [MPM,M]=ptTrackLinker(M)
%
% INPUT         M          : M stack as returned by the tracker functions
%                                  M = [y x y x]   [y x y x]   [y x y x]
%                                         ...    ,    ...    ,    ...
%                                       t1   t2     t2   t3     t3   t4
%                                          1           2           3  
% OUTPUT        MPM        : Magic Position Matrix 
%                                MPM = [ y  x  y  x  y  x ... ]
%                                         t1    t2    t3
%               M          : Rearranged M matrix.
%
% DEPENDENCES   fsmTrackLinker uses { }
%               fsmTrackLinker is used by { ptTrackCells } 
%
% Aaron Ponti, 2002
%
% altered from fsmTrackLinker by Colin Glass, Feb 2004
%
% the changes in this code result in the fact, that MPM returns only one
% tracked thing per row. Original code can return more than one per row,
% seperated by zeros

% Let the user know we're starting to link
fprintf (1, 'Starting track linkage process...\n');

% Initialize counter and waitbar
counter = 0;
tot = 2 * (size(M,3)-1);
%h=waitbar(0,'Linking...');

% Reorganize M
for counter1 = 1 : size(M,3) - 1
  
   % counter
   counter = counter + 1;
   fprintf(1, 'Linking frame # %d ...\n', counter);
   
   % Read speckle positions at time point (=img) counter1
   start = (M(:, 3:4, counter1));
   stop = (M(:, 1:2, counter1+1));
   
   % Re-arrange stop (and therefore M) to correspond to the sequence of start
   tM = zeros(size(start, 1), 4);
   
   for counter2 = 1:size(stop, 1)

      %Attention: this if /else statement is the only difference to Aarons original code 
      if start(counter2, 1) ~= 0 | start(counter2, 2) ~= 0 
                        
         t = start(counter2, 1) == stop(:,1);
         u = start(counter2, 2) == stop(:,2);
         y = find(t & u);
                        
         %
         % ANALYSIS
         %
                        
         % No matching found -> error!
         if isempty(y)
             fprintf(1, 'Time points %d to %d.\n', counter1, counter1 + 1);
             warning('ptTrackLinker: Warning! Correspondance not found.');
             tM(counter2,:)=0; % -1;   
%               if M(counter2,3,counter1+1)~=0
%                        M(end+1,:,:)=0;
%                        tM(end+1,:)=0;
%                        tM(end,3:4)=M(counter2,3:4,counter1+1);
%              end
         end
                        
         % Only one entry found in stop
         if length(y)==1
             tM(counter2,:)=M(y,:,counter1+1);
             stop(y,:)=-3;
         end
                        
         % More than one entry found, but either 'no speckle' (0) 
         %   or 'already treated' (-3)
         if length(y)>1 & (start(counter2,1)~=0 | start(counter2,1)~=-3)
             tM(counter2,:)=M(y(1),:,counter1+1);
             stop(y(1),:)=-3;
         end
                        
         % More than one repetition of a speckle found (~=0 & ~=-3)
         if  length(y)>1 & start(counter2,1)~=0 & start(counter2,1)~=-3
             error('fsmTrackLinker: Warning! Not all repetitions have been removed.');
         end        
         
         % 
         % END OF ANALYSIS
         %
      else
         tM(counter2,:)=0;
                
         if (M(counter2,3,counter1+1)~=0 | M(counter2,4,counter1+1)~=0) & (M(counter2,1,counter1+1)==0 & M(counter2,2,counter1+1)==0)
            M(end+1,:,:)=0;
            %=M(end,3:4,counter1);
            tM(end+1,:)=0;
            tM(end,3:4)=M(counter2,3:4,counter1+1);
         end
      end
   end
    
   % Replace M with re-ordered one
   M(:,:,counter1+1)=tM;
    
   % Reset tM
   tM=zeros(size(tM));
    
   % Update wait bar
   %waitbar(counter/tot,h);   
   %fprintf(1, 'Linking # frame %d ...\n', counter);
end

MPM(:,1:2)=M(:,1:2,1);
counter = 0;

% Keep the user up-to-date on what is happening
fprintf (1, 'Postprocessing linked frames...\n');

% Remove info that is not needed anymore
for counter3=2:size(M,3)
    
   % counter
   %counter = counter + 1;
    
   MPM(:,(counter3-1)*2+(1:2))=M(:,1:2,counter3);
    
   % Update wait bar
   %waitbar(counter/tot,h);
   %fprintf(1, 'Linking frame # %d ...\n', counter);
end

MPM(:,counter3*2+(1:2))=M(:,3:4,counter3);

% Let the user know we've finished
fprintf (1, 'Finished linking tracks!\n');
% Close waitbar
%close(h);
