function stuffplotter(hObject)

handles=guidata(hObject);


start=round((handles.postpro.analfirstimg- handles.jobvalues.firstimage)/handles.jobvalues.increment)+1;

stop = floor((handles.postpro.anallastimg - handles.jobvalues.firstimage)/handles.jobvalues.increment+0.00001)+1 ;



 PROPPERIMG=zeros(6,start-stop+1);



for pic=start:stop
    
    whatcells=[];
    
	if ~isempty(handles.selectedcells)
        whatcells=zeros(size(handles.selectedcells,1),2);
	else
        whatcells=zeros(size(handles.MPM,1),2);
	end

    
    
    if ~isempty(handles.selectedcells) 
        whatcells(:,:)=handles.MPM(handles.selectedcells,2*pic-1:2*pic);
    else  
        whatcells(:,:)=handles.MPM(:,2*pic-1:2*pic);
    end
    
   whatcells=whatcells(find(whatcells(:,1)&whatcells(:,2)),:); 
   
   PROPPERIMG(6,pic)=size(whatcells,1);
    
    properRows= ismember(round(handles.PROPERTIES(:,1:2,pic)),round(whatcells),'rows');
    
    takeintoaccount=handles.PROPERTIES(properRows,:,pic);
 
    
    %how many clusters
    
    [uniqclust,rows,invert]=unique(takeintoaccount(:,3));
    PROPPERIMG(1,pic)=length(uniqclust);
    
    
    %how many cells per cluster
    PROPPERIMG(2,pic)=sum(takeintoaccount(rows,4))/PROPPERIMG(1,pic);
    
    %how big an area per cell
    
    PROPPERIMG(3,pic)=sum(takeintoaccount(:,5))/PROPPERIMG(6,pic);
    
    %how big an area per cluster
    bridge=takeintoaccount(:,5).*(takeintoaccount(:,4));
    PROPPERIMG(4,pic)=sum(bridge(rows))/length(rows);
    
    %perimeter/area
    uff=find(takeintoaccount(:,6));
    PROPPERIMG(5,pic)=sum(takeintoaccount(uff,6))/length(uff);
    
   
end 

if get(handles.GUI_ad_numberofthings_rb,'Value')
    
        figure
        
        subplot(2,2,1);plot(PROPPERIMG(1,:)), title('how many clusters')
        xlabel('Frames')
        ylabel('# of clusters')

        subplot(2,2,2);plot(PROPPERIMG(2,:)), title('how many cells per cluster')
        xlabel('Frames')
        ylabel('cells per cluster')
        
        subplot(2,2,3);plot(PROPPERIMG(6,:)), title('how many cells')
        xlabel('Frames')
        ylabel('# of cell')
        
end   
        
if get(handles.GUI_ad_areas_rb,'Value')
        
        figure
        
        subplot(1,2,1);plot(PROPPERIMG(3,:)), title('how big an area per cell')
        xlabel('Frames')
        ylabel('area per cell')
        
        subplot(1,2,2); plot(PROPPERIMG(4,:)), title('how big an area per cluster')
        xlabel('Frames')
        ylabel('area per cluster')
end

if get(handles.GUI_ad_perimeter_rb,'Value')
        
        figure, plot(PROPPERIMG(5,:)), title('perimeter/area of clusters')
        xlabel('Frames')
        ylabel('perimeter/area')
end
    
%     
%       figure, plot(numbercells),title('number of cells')
%     xlabel('Frames')
%     ylabel('# cells')
%     
%     
%     figure, 
%     subplot(2,1,1); plot(avardist),title('avarage distance travelled by a cell')
%     xlabel('Frames')
%     ylabel('distance per frame (in pixel)')
%   
%     subplot(2,1,2); plot(speedvar),title('variance of avarage distance travelled by a cell')
%       xlabel('Frames')
%     ylabel('distance per frame (in pixel)')