function movieMaker(hObject) 

handles=guidata(hObject);



[numrows,numcols]=size(handles.MPM);

ImageNamesList = handles.jobvalues.imagenameslist;

dragTailLength = handles.postpro.dragtail;

start = round((handles.postpro.moviefirstimg- handles.jobvalues.firstimage)/handles.jobvalues.increment)+1
if start < dragTailLength+2
    start = dragTailLength+2;
end

stop = floor((handles.postpro.movielastimg - handles.jobvalues.firstimage)/handles.jobvalues.increment+0.00001)+1 

if stop > handles.jobvalues.lastimage
    stop = handles.jobvalues.lastimage;
end

imagedirectory = handles.jobvalues.imagedirectory;
saveallpath = handles.postpro.saveallpath;

cd(saveallpath)

makeQTmovie('start','trackmov.mov')


for movieStep = start:stop
    
    whatcells=[];
    
	if ~isempty(handles.selectedcells)
        whatcells=zeros(size(handles.selectedcells,1),2);
           whatcells(:,:)=handles.MPM(handles.selectedcells,(2*movieStep-1):(2*movieStep));
        
	else
        whatcells=zeros(size(handles.MPM,1),2);
        whatcells(:,:)=handles.MPM(:,(2*movieStep-1):(2*movieStep));
	end

    
    
    %whatcells defines which cells will be taken into account. This is
    %either defined by the user (PolyTrack_PP) or it will just be all cells
    %within the current picture
     
     whatcells=find(whatcells(:,1)&whatcells(:,2)); 
    
     cd(imagedirectory)

     
     name = char(ImageNamesList(movieStep)); 
     nowImg=imreadnd2(name,0,handles.jobvalues.intensityMax);

		
		
		
		
		[rows,cols]= find(handles.MPM);
		
		
		
		if ~isempty(handles.postpro.figureSize)
		     figure('Position',handles.postpro.figureSize) 
             imshow(nowImg,[])
        else
             figure, imshow(nowImg,[])
        end
        
		hold on
		
		
		%one colour per time step
		cmap=jet(dragTailLength+1);
		%cmap=jet(numcols/2-1);
		counter=0;
		for i=(2*(movieStep-dragTailLength)):2:(2*movieStep)
            
                counter=counter+1;
                vec=handles.MPM(whatcells,i-3:i);
                [rows,cols]=find(vec==0);
                rows=unique(rows);
                vec(rows,:)=0;
                
                
                 for h=1:size(vec,1)
                     
                 
                     %   if vec(h,1)~=0 & sqrt((vec(:,1)-vec(:,3)).^2 + (vec(:,2)-vec(:,4)).^2)<80
                             %   plot(vec(h,1),vec(h,2),'.')
                             ph=[];
                              ph=  plot(vec(h,1:2:3),vec(h,2:2:4));
                              set(ph,'Color',cmap(counter,:));
                              clear ph
                              %  end
                  
                end
                
                
                
		end
		plot(vec(:,3),vec(:,4),'r.')
		
		hold off
		
        
        cd(saveallpath)
  
		makeQtmovie('addaxes',gca);
		close
end

makeQtmovie('finish');











%         
%         
% % % % %         for h=1:size(vec,1)
% % % % %          
% % % % %                 if vec(h,1)~=0
% % % % %                      %   plot(vec(h,1),vec(h,2),'.')
% % % % %                      ph=[];
% % % % %                       ph=  plot(vec(h,1:2:3),vec(h,2:2:4));
% % % % %                       set(ph,'Color',cmap(h,:));
% % % % %                       clear ph
% % % % %                 end
% % % % %           
% % % % %         end








% % hold on
% % 
% % PROPPERIMG=getinfor(PROPERTIES)
% % 
% % figure,title('coordinates per image')
% % hold on
% % 
% % PROPERTIES

% 	cmap=jet(numcols/2-1);
% 	
% 	counter=1
% 	
% 	w = waitforbuttonpress;
% 	while w==0 & (counter < numrows+1)
%            f=figure, imshow(pics(end).info)
%            
%             onecell=handles.MPM(counter,:);
%             index=find(onecell);
%              beginn=min(index);
%              nomore=max(index);
%         
%             if nomore-beginn >5
%                 
%                  
%              for h=1:2:((nomore-beginn)-1)
%                        hold on
%                          %   plot(vec(h,1),vec(h,2),'.')
%                          ph=[];
%                           ph=  plot(onecell(beginn+h-1:2:beginn+h+1),onecell(beginn+h:2:beginn+h+2));
%                           set(ph,'Color',cmap((beginn+h)/2,:));
%                           clear ph
%                            hold off
%               
%             end
            
            
            
% %             
% %             hold on
% %             handl(counter)=plot(onecell(beginn:2:nomore-1),onecell(beginn+1:2:nomore));
% %            set(handl(counter),'Color',cmap(counter,:));
% %            hold off
%          end
%          
%        w = waitforbuttonpress;
%        
%       close
% 
%       counter=counter+1    
%    end

% 
% %make a nice little picture with the paths of the cells in it
% 
% %one colour per cell
% 
% ;
% for p=1:numrows
%     
%     onecell=handles.MPM(p,:);
%     index=find(onecell);
%     beginn=min(index);
%     nomore=max(index);
%     
%     if nomore-beginn >3
%             handl(p)=plot(onecell(beginn:2:nomore-1),onecell(beginn+1:2:nomore));
%            set(handl(p),'Color',cmap(p,:));
%    end
%    
% end 
% 
% 

