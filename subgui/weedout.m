
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



[firsts,cols]= find(handles.MPM);
firsts=unique(firsts);
handles.MPM=handles.MPM(firsts,:);





%The goal of this little routine is to relink tracks
for i=1:size(handles.MPM,1)
        onecell=handles.MPM(i,:);
        index=find(onecell);
        beginn(i)=min(index);
        nomore(i)=max(index);
end

counter=0;
while counter<size(handles.MPM,1)-0.3
       
    counter=counter+1;
    
    
    near= find((nomore(counter)-beginn(counter+1:end))<minusframes & (nomore(counter)-beginn(counter+1:end))>0 & (nomore(counter)-nomore(counter+1:end))<plusframes);
    
    if ~isempty (near)
        if nomore(counter)+1<size(handles.MPM,2)
               [distance,chuck] = min(sqrt((handles.MPM(counter,nomore(counter)-1)-handles.MPM(near(:)+...
                        counter,nomore(counter)+1)).^2+(handles.MPM(counter,nomore(counter))...
                        -handles.MPM(near(:)+counter,nomore(counter)+2)).^2));
             
                if distance < maxdistpostpro
                    
                    handles.MPM(counter,nomore(counter)+1:end)=handles.MPM(chuck+counter,nomore(counter)+1:end);
                    
                    handles.MPM(chuck+counter,:)=0;
                    
                end
       end

    end

end



clear chuck;

clear firsts;
[firsts,cols]= find(handles.MPM);
firsts=unique(firsts);
handles.MPM=handles.MPM(firsts,:);


[rows,cols]=find(handles.MPM);





% here we chuck out cells that didn't stay longer than mintrack frames
sorrows = sort(rows);
[uniqueEntries, uniqueIdx] = unique(sorrows);
%uniqueIdx returns the last occurence of the respective unique entry
%having sorted rows before, we can now count the number of occurences
if size(uniqueEntries,1) > size(uniqueEntries,2);
        uniqueIdx = [0;uniqueIdx];
else
        uniqueIdx = [0,uniqueIdx];
end 
numberOfOccurences = diff(uniqueIdx); 
chuck = uniqueEntries(find(numberOfOccurences < minimaltrack*2));
handles.MPM(chuck,:) = [];

%set(handles.GUI_app_autopostpro_cb,'Value',1);

guidata(hObject,handles);

MPM=handles.MPM
saveallpath=get(handles.GUI_fm_saveallpath_ed,'String')

cd(saveallpath)
save('MPM', 'MPM')