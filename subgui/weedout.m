
function weedout(hObject)

handles=guidata(hObject);

% 
% if get(handles.GUI_app_autopostpro_cb,'Value')==1
%    set(handles.GUI_app_autopostpro_cb,'Value',0);
% end
%     
    
minimaltrack= handles.postpro.minimaltrack;
maxdistpostpro= handles.postpro.maxdistpostpro ;

minusframes = handles.postpro.minusframes;
plusframes = handles.postpro.plusframes;


%if there are totally empty rows, erase them
[firsts,cols]= find(handles.MPM);
firsts=unique(firsts);
handles.MPM=handles.MPM(firsts,:);





%The goal of this little routine is to relink tracks
%first we find (for every cell) the frame in which it appears and the frame
%in which it disappears
for i=1:size(handles.MPM,1)
        onecell=handles.MPM(i,:);
        index=find(onecell);
        beginn(i)=min(index);
        nomore(i)=max(index);
end

counter=0;
%loop over the number of cells. Note: the older the cell is, the lower 
%it's row index. So we always try to add information of cells with higher
%row index (newer cells) to ones with lower row index (older cells) and 
%then erase the newer cell. 
%in this way we reduce the likelyhood of asigning a new number to a cell
%already found 
while counter < (size(handles.MPM,1)-0.3)
       
    linked = 0;
    
    counter=counter+1;
    
    
    %we look for tracks which begin only a few frames before our current
    %track stops and try to link 
    near= find(  (nomore(counter)-beginn(counter+1:end)) < minusframes ...
               & (nomore(counter)-beginn(counter+1:end)) > 0 ...
               & (nomore(counter)-nomore(counter+1:end)) < plusframes ...
               &  beginn(counter+1:end) > 1.1 ...
               & (nomore(counter)+1.9) < size(handles.MPM,2)  );
    
    if ~isempty (near)
    
           [distance,chuck] = min( sqrt( ( handles.MPM(counter , nomore(counter)-1) - handles.MPM(near(:)+counter , nomore(counter)+1) ).^2 ... 
                                           + ( handles.MPM(counter , nomore(counter)) - handles.MPM(near(:)+counter , nomore(counter)+2) ).^2  ) );
         
            if distance < maxdistpostpro
                
                %add the found track to the current track
                handles.MPM(counter , (nomore(counter)+1):end) = handles.MPM(near(chuck)+counter , (nomore(counter)+1):end);
                
                %erase the redundant track
                handles.MPM(near(chuck)+counter , :) = 0;
                
                %make sure the current track gets processed again, with the new stop location. Maybe we can find a new link there 
                nomore(counter) = nomore(near(chuck)+counter);
                
                %erase begin and stop location of our allocated track (now
                %part of the current track and no longer an individual
                %track)
                nomore(near(chuck)+counter) = 0;
                beginn(near(chuck)+counter) = 0;
                
                counter = counter - 1;
                
                linked = 1;
                
            end
    
    end
    
    distance=[];
    chuck=[];
    
    %we look for tracks which begin two frames after our current
    %track stops and try to link
    if ~linked
          gaps = find(  (nomore(counter)+3.9) < size(handles.MPM,1) ...
                      & ((nomore(counter)-beginn(counter+1:end)) == -2) )
                    
          if ~isempty (gaps)
    
                [distance,chuck] = min( sqrt( ( handles.MPM(counter , nomore(counter)-1) - handles.MPM(gaps(:)+counter , nomore(counter)+3) ).^2 ... 
                                           + ( handles.MPM(counter , nomore(counter)) - handles.MPM(gaps(:)+counter , nomore(counter)+4) ).^2  ) );
                
                if distance < maxdistpostpro
                        
                        %add the found track to the current track
                        handles.MPM(counter , (nomore(counter)+3):end) = handles.MPM(gaps(chuck)+counter , (nomore(counter)+3):end);
                        
                        %fill in the gap
                        handles.MPM(counter , nomore(counter)+1)=   handles.MPM(counter , nomore(counter)-1) +...
                                                                    round( (handles.MPM(counter , nomore(counter)+3) - handles.MPM(counter , nomore(counter)-1))/2 );
                        handles.MPM(counter , nomore(counter)+2)=   handles.MPM(counter , nomore(counter)) +...
                                                                    round( (handles.MPM(counter , nomore(counter)+4) - handles.MPM(counter , nomore(counter)))/2 );
                        %erase the redundant track
                        handles.MPM(gaps(chuck)+counter , :) = 0;
                        
                        %make sure the current track gets processed again, with the new stop location. Maybe we can finf a new link there 
                        nomore(counter) = nomore(gaps(chuck)+counter);
                        
                        %erase begin and stop location of our allocated track (now
                        %part of the current track and no longer an individual
                        %track)
                        nomore(gaps(chuck)+counter) = 0;
                        beginn(gaps(chuck)+counter) = 0;
                        
                        counter = counter - 1;
                end  
                        
          end
        
     end
        
end


clear chuck;

%again erase all tracks with no entries at all
clear firsts;
[firsts,cols]= find(handles.MPM);
firsts=unique(firsts);
handles.MPM=handles.MPM(firsts,:);


% here we chuck out cells that didn't stay longer than mintrack frames
[rows,cols]=find(handles.MPM);

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
handles.MPM(chuck,:) = [];

%set(handles.GUI_app_autopostpro_cb,'Value',1);


%save the altered data
guidata(hObject,handles);

MPM=handles.MPM
saveallpath=get(handles.GUI_fm_saveallpath_ed,'String')

cd(saveallpath)
save('MPM', 'MPM')