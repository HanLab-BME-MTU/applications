
function MPM = weedout(MPM,plusframes,minusframes,maxdistpostpro,minimaltrack,saveallpath)
% weedout filters and relinks the tracks given in MPM
%
% SYNOPSIS       MPM = weedout(MPM,plusframes,minusframes,maxdistpostpro,minimaltrack,saveallpath)
%
% INPUT          MPM        : Magic Position Matrix 
%                               MPM = [ y  x  y  x  y  x ... ]
%                                        t1    t2    t3
%                plusframes : how many frames after the currently analysed
%                             track stops must a candidate for linkage survive
% 		  		 minusframes : how many frames (at most!) before the currently analysed
%                              track stops can a candidate for linkage be
%                              present
% 				 maxdistpostpro : max distance for relinking
% 				 minimaltrack : minimal length of tracks (ales they get
% 			                	erased
% 				 saveallpath : where shall the new MPM be saved
%
%
% OUTPUT         altered MPM          
%
% DEPENDENCIES   weedout  uses {nothing}
%                                  
%                weedout is used by { PolyTrack_PP }
%
% Colin Glass, Feb 04



% % % 
% % % if get(handles.GUI_app_autopostpro_cb,'Value')==1
% % %    set(handles.GUI_app_autopostpro_cb,'Value',0);
% % % end
% % %     

%if there are totally empty rows, erase them
[firsts,cols]= find(MPM);
firsts=unique(firsts);
MPM=MPM(firsts,:);





%The goal of this little routine is to relink tracks
%first we find (for every cell) the frame in which it appears and the frame
%in which it disappears
for i=1:size(MPM,1)
        onecell=MPM(i,:);
        index=find(onecell);
        appear(i)=min(index);
        disappear(i)=max(index);
end

counter=0;
%loop over the number of cells. Note: the older the cell is, the lower 
%it's row index. So we always try to add information of cells with higher
%row index (newer cells) to ones with lower row index (older cells) and 
%then erase the newer cell. 
%in this way we reduce the likelyhood of asigning a new number to a cell
%already found 
while counter < (size(MPM,1)-0.3)
       
    linked = 0;
    
    counter=counter+1;
    
    
    %we look for tracks which begin only a few frames before our current
    %track stops and try to link 
    near= find(  (disappear(counter)-appear(counter+1:end)) < minusframes ...
               & (disappear(counter)-appear(counter+1:end)) > 0 ...
               & (disappear(counter)-disappear(counter+1:end)) < plusframes ...
               &  appear(counter+1:end) > 1.1 ...
               & (disappear(counter)+1.9) < size(MPM,2)  );
    
    if ~isempty (near)
    
           [distance,chuck] = min( sqrt( ( MPM(counter , disappear(counter)-1) - MPM(near(:)+counter , disappear(counter)+1) ).^2 ... 
                                           + ( MPM(counter , disappear(counter)) - MPM(near(:)+counter , disappear(counter)+2) ).^2  ) );
         
            if distance < maxdistpostpro
                
                %add the found track to the current track
                MPM(counter , (disappear(counter)+1):end) = MPM(near(chuck)+counter , (disappear(counter)+1):end);
                
                %erase the redundant track
                MPM(near(chuck)+counter , :) = 0;
                
                %make sure the current track gets processed again, with the new stop location. Maybe we can find a new link there 
                disappear(counter) = disappear(near(chuck)+counter);
                
                %erase begin and stop location of our allocated track (now
                %part of the current track and no longer an individual
                %track)
                disappear(near(chuck)+counter) = 0;
                appear(near(chuck)+counter) = 0;
                
                counter = counter - 1;
                
                linked = 1;
                
            end
    
    end
    
    distance=[];
    chuck=[];
    
    %we look for tracks which begin two frames after our current
    %track stops and try to link
    if ~linked
          gaps = find(  (disappear(counter)+3.9) < size(MPM,1) ...
                      & ((disappear(counter)-appear(counter+1:end)) == -2) )
                    
          if ~isempty (gaps)
    
                [distance,chuck] = min( sqrt( ( MPM(counter , disappear(counter)-1) - MPM(gaps(:)+counter , disappear(counter)+3) ).^2 ... 
                                           + ( MPM(counter , disappear(counter)) - MPM(gaps(:)+counter , disappear(counter)+4) ).^2  ) );
                
                if distance < maxdistpostpro
                        
                        %add the found track to the current track
                        MPM(counter , (disappear(counter)+3):end) = MPM(gaps(chuck)+counter , (disappear(counter)+3):end);
                        
                        %fill in the gap
                        MPM(counter , disappear(counter)+1)=   MPM(counter , disappear(counter)-1) +...
                                                                    round( (MPM(counter , disappear(counter)+3) - MPM(counter , disappear(counter)-1))/2 );
                        MPM(counter , disappear(counter)+2)=   MPM(counter , disappear(counter)) +...
                                                                    round( (MPM(counter , disappear(counter)+4) - MPM(counter , disappear(counter)))/2 );
                        %erase the redundant track
                        MPM(gaps(chuck)+counter , :) = 0;
                        
                        %make sure the current track gets processed again, with the new stop location. Maybe we can finf a new link there 
                        disappear(counter) = disappear(gaps(chuck)+counter);
                        
                        %erase begin and stop location of our allocated track (now
                        %part of the current track and no longer an individual
                        %track)
                        disappear(gaps(chuck)+counter) = 0;
                        appear(gaps(chuck)+counter) = 0;
                        
                        counter = counter - 1;
                end  
                        
          end
        
     end
        
end


clear chuck;

%again erase all tracks with no entries at all
clear firsts;
[firsts,cols]= find(MPM);
firsts=unique(firsts);
MPM=MPM(firsts,:);


% here we chuck out cells that didn't stay longer than mintrack frames
[rows,cols]=find(MPM);

sorrows = sort(rows);
[uniqueEntries, uniqueIdx] = unique(sorrows);
%uniqueIdx returns the last occurence of the respective unique entry
%having sorted rows before, we can now count the number of occurences\

if size(uniqueEntries,1) > size(uniqueEntries,2);
        uniqueIdx = [0;uniqueIdx];
else
        uniqueIdx = [0,uniqueIdx];
end 

numberOfOccurences = diff(uniqueIdx); 
chuck = uniqueEntries(find(numberOfOccurences < minimaltrack*2-1));
MPM(chuck,:) = [];

%set(handles.GUI_app_autopostpro_cb,'Value',1);




saveallpath=get(handles.GUI_fm_saveallpath_ed,'String')

cd(saveallpath)
save('MPM', 'MPM')