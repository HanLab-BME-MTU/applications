function speed(hObject)

handles=guidata(hObject);


[numrows,numcols]=size(handles.MPM)

dist=zeros(numrows,1);

for i=1:2:numcols-3
        
        vec=handles.MPM(:,i:i+3);
        [rows,cols]=find(vec==0);
        rows=unique(rows);
        vec(rows,:)=0;
        
        numbercells((i+1)/2,1)=(numrows-length(rows))
        dist(:,(i+1)/2)=sqrt((vec(:,1)-vec(:,3)).^2 + (vec(:,2)-vec(:,4)).^2);
        avardist((i+1)/2,1)=sum(dist(1:numrows,(i+1)/2))/(numrows-length(rows));
        
        realdist=dist(find(dist(:,(i+1)/2)),(i+1)/2);
        speedvar((i+1)/2)=var(realdist);
         
         
    end
    
    
    
    figure, plot(numbercells),title('number of cells')
    xlabel('Frames')
    ylabel('# cells')
    
    
    figure, 
    subplot(2,1,1); plot(avardist),title('avarage distance travelled by cells')
    xlabel('Frames')
    ylabel('distance per frame (in pixel)')
  
    subplot(2,1,2); plot(speedvar),title('variance of avarage distance travelled by a cell')
      xlabel('Frames')
    ylabel('variance')
    
   
        
saveallpath=get(handles.GUI_fm_saveallpath_ed,'String')

cd(saveallpath)
save('distances', 'dist')
save('speedvariances','speedvar')
save('avardist','avardist')
    
   guidata(hObject, handles);