function ptPlotCellValues.m(hObject)
% ptPlotCellValues.m plots information gathered in cellProps. Certain details
%              get derived from cellProps first
%
% SYNOPSIS       ptPlotCellValues.m(hObject)
%
% INPUT          hObject : handle to an object within PolyTrack_PP
%                
% OUTPUT         altered MPM          
%
% DEPENDENCIES   ptPlotCellValues.m  uses {nothing}
%                                  
%                ptPlotCellValues.m is used by { PolyTrack_PP }
%
% Colin Glass, Feb 04


%We now want to rearrange certain data (cellProps), so that we can plot it.  

handles=guidata(hObject);


start = round((handles.postpro.analfirstimg- handles.jobvalues.firstimage)/handles.jobvalues.increment)+1;
if start < 1
    start = 1;
end

stop = floor((handles.postpro.anallastimg - handles.jobvalues.firstimage)/handles.jobvalues.increment+0.00001)+1;

if stop > handles.jobvalues.lastimage
    stop = handles.jobvalues.lastimage;
end

saveallpath=handles.postpro.saveallpath;


NewProps=zeros(6,start-stop+1);


%erase this information. it was calculated before we cleaned up the results
handles.cellProps(:,4,:)=0;
cellProps = handles.cellProps;


  whatcells=[];
    
	if ~isempty(handles.selectedcells)
        whatcells = zeros(size(handles.selectedcells,1),2);
        handles.singlecells = zeros(size(handles.selectedcells,1),1)
        handles.clustercells = zeros(size(handles.selectedcells,1),1)
        
	else
        whatcells = zeros(size(handles.MPM,1),2);
        handles.singlecells = zeros(size(handles.MPM,1),2);
        handles.clustercells = zeros(size(handles.MPM,1),2);
        
	end

    



for pic=start:stop
 
    %now we quickly calculate the number of cells within a cluster. This
    %may seem a little bit late to do so, but only now we know which cells
    %are actually accepted by postprocessing
    
	belo = sort(handles.cellProps(:,3,pic));
	[uniqueEntries, uniqueIdx] = unique(belo);
	%uniqueIdx returns the last occurence of the respective unique entry
	%having sorted m before, we can now count the number of occurences
	if size(uniqueEntries,1) > size(uniqueEntries,2);
            uniqueIdx = [0;uniqueIdx];
	else
            uniqueIdx = [0,uniqueIdx];
	end 
	numberOfOccurences = diff(uniqueIdx); 

    %calculate the number of members for every cluster
     for indUniEnt = 1:length(uniqueEntries)
           order = ismember(handles.cellProps(:,3,pic),uniqueEntries(indUniEnt));
           handles.cellProps(order ,4,pic) = numberOfOccurences(indUniEnt);
     end
     %done
     
     
     
    whatcells=[];
    %if the user has selected certain cells, only use these
	if ~isempty(handles.selectedcells)
        whatcells = zeros(size(handles.selectedcells,1),2);
      
	else
        whatcells = zeros(size(handles.MPM,1),2);
       
	end

    
    if ~isempty(handles.selectedcells) 
        whatcells(:,:)=handles.MPM(handles.selectedcells,(2*pic-1):(2*pic));
       % presentCells(:,:)=handles.MPM(:,(2*pic-1):(2*pic));
        
    else  
        whatcells(:,:)=handles.MPM(:,(2*pic-1):(2*pic));
       % presentCells(:,:)=handles.MPM(:,(2*pic-1):(2*pic));
    end
    
    %whatcells defines which cells will be taken into account. This is
    %either defined by the user (PolyTrack_PP) or it will just be all cells
    %within the current picture
    good = find(whatcells(:,1) | whatcells(:,2));
    whatcells = whatcells(good,:); 

    
%     presentCells = presentCells(find(presentCells(:,1) & presentCells(:,2)),:); 
%     properRows = ismember(round(handles.cellProps(:,1:2,pic)),round(presentCells),'rows');
%     referenceAllCells = handles.cellProps(properRows,:,pic);
%     
%     clear properRows
    


    %here we find all the rows of cellProps that correspond to cells we are
    %intersted in. cellProps(:,1:2,pic) and whatcells are both [x,y]
    %coordinates
    
    [properRows,indWhatCells] = ismember(round(handles.cellProps(:,1:2,pic)),round(whatcells),'rows');
    NewProps(6,pic) = size(find(properRows),1);
    takeIntoAccount = handles.cellProps(properRows,:,pic);
    
    indWhatCells = indWhatCells(find(indWhatCells));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% REMINDER:
	% one row of cellProps gives information on one set of coordinates (one
	% cell). takeIntoAccount has the same layout
	% cellProps=zeros(length(coord),6);
	% cellProps(:,1)=coord(:,1);
	% cellProps(:,2)=coord(:,2);
	% cellProps(:,3)=belongsto(:);  (number of cluster - label)
	% cellProps(:,4)=numberOfOccurences(:);  (how many cells in the cluster
	%                                          this cell is in)
	% cellProps(:,5)=bodycount(:);  (area of the cluster with the number given in belongsto)
	% cellProps(:,6)=perimdivare(:);  (cluster)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    
    %how many clusters. At least two members.
    Clusters = find(takeIntoAccount(:,4)>1.5);
    
     handles.clustercells(1:length(Clusters),pic) = indWhatCells(Clusters);
    
    [uniqclust,uniCluRows,invert]=unique(takeIntoAccount(Clusters,3));
    NewProps(1,pic)=length(uniqclust);
    
    
    %how many cells per cluster
    NewProps(2,pic)=sum(takeIntoAccount(Clusters(uniCluRows),4))/NewProps(1,pic);
    
    %how big an area per cell
    bridge=takeIntoAccount(:,5)./takeIntoAccount(:,4);
    NewProps(3,pic)=sum(bridge)/NewProps(6,pic);
    
    %how big an area per cluster
  
    NewProps(4,pic)=sum(takeIntoAccount(Clusters(uniCluRows),5))/length(uniCluRows);
    
    %perimeter/area
    NewProps(5,pic)=sum(takeIntoAccount(Clusters(uniCluRows),6))/length(uniCluRows);
    
    %percentage single cells
    indSingle = find(takeIntoAccount(:,4)==1);
    
    handles.singlecells (1:length(indSingle),pic) = indWhatCells(indSingle);
    
    NewProps(7,pic) = size(indSingle,1) / NewProps(6,pic);
    
   
end 

% if get(handles.GUI_ad_numberofthings_rb,'Value')
%     
%         h_fig=figure,title(handles.jobvalues.imagename);
%         
%         ymax=max(NewProps(1,:));
%         subplot(2,2,1);plot(NewProps(1,:)), title('how many clusters');
%         xlabel('Frames');
%         ylabel('# of clusters');
%         axis([0 stop 0 ymax]);
%         
%         ymax=max(NewProps(2,:));
%         subplot(2,2,2);plot(NewProps(2,:)), title('how many cells per cluster');
%         xlabel('Frames');
%         ylabel('cells per cluster');
%         axis([0 stop 0 ymax]);
%  
%         ymax=max(NewProps(6,:));
%         subplot(2,2,3);plot(NewProps(6,:)), title('how many cells');
%         xlabel('Frames');
%         ylabel('# of cell');
%         axis([0 stop 0 ymax]);
% 
%         
%         
% 
% 
% 		
% 		hgsave(h_fig,[saveallpath filesep 'numberCellsClust.fig']);
% 		print(h_fig, [saveallpath filesep 'numberCellsClust.eps'],'-depsc2','-tiff');
% 		print(h_fig, [saveallpath filesep 'numberCellsClust.tif'],'-dtiff');
% 		
% 
%         
% end   
%         
% if get(handles.GUI_ad_areas_rb,'Value')
%         
%        h_fig=figure,title(handles.jobvalues.imagename);
%         
%        ymax=max(NewProps(3,:));
%         subplot(1,2,1);plot(NewProps(3,:)), title('how big an area per cell');
%         xlabel('Frames');
%         ylabel('area per cell');
%         axis([0 stop 0 ymax]);
%         
%         ymax=max(NewProps(4,:));
%         subplot(1,2,2); plot(NewProps(4,:)), title('how big an area per cluster');
%         xlabel('Frames');
%         ylabel('area per cluster');
%         axis([0 stop 0 ymax]);
%         
%         
% 		hgsave(h_fig,[saveallpath filesep 'areas.fig']);
% 		print(h_fig, [saveallpath filesep 'areas.eps'],'-depsc2','-tiff');
% 		print(h_fig, [saveallpath filesep 'areas.tif'],'-dtiff');
% 		
%         
% end

% if get(handles.GUI_ad_perimeter_rb,'Value')
%         
%         h_fig=figure,title(handles.jobvalues.imagename);
%         ymax=max(NewProps(5,:));
%         subplot(1,1,1); plot(NewProps(5,:)), title('perimeter/area of clusters');
%         xlabel('Frames');
%         ylabel('perimeter/area');
%         axis([0 stop 0 ymax]);
%         
%            
% 		hgsave(h_fig,[saveallpath filesep 'areasPerim.fig']);
% 		print(h_fig, [saveallpath filesep 'areasPerim.eps'],'-depsc2','-tiff');
% 		print(h_fig, [saveallpath filesep 'areasPerim.tif'],'-dtiff');
% 		
%         
%         
% end
%     

percentageSingle = NewProps(7,:);

% h_fig=figure,title(handles.jobvalues.imagename);
% ymax=max(percentageSingle);
% subplot(1,1,1); plot(percentageSingle), title('Percentage of single cells (red - filtered, blue - non-filtered)');
% hold on;
% plot(percentageSingleFiltered, 'r');
% xlabel('Frames');
% ylabel('percentage of single cells');
% axis([0 stop 0 1]);
% hold off;
%         
%            
% hgsave(h_fig,[saveallpath filesep 'percentSingle.fig']);
% print(h_fig, [saveallpath filesep 'percentSingle.eps'],'-depsc2','-tiff');
% print(h_fig, [saveallpath filesep 'percentSingle.tif'],'-dtiff');
	
        
guidata(hObject, handles);

cd(saveallpath);
save ('NewProps', 'NewProps');
save ('percentageSingle','percentageSingle');
