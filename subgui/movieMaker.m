function movieMaker(hObject) 
% movieMaker makes movies from the information gathered in the amin
%            analysis
%
% SYNOPSIS       movieMaker(hObject)
%
% INPUT          hObject : handle to an object within PolyTrack_PP
%                MPM
%                ImageNamesList
%                dragTailLength...
%
% OUTPUT         saves movies to disk          
%
% DEPENDENCIES   movieMaker  uses {nothing}
%                                  
%                movieMaker is used by { PolyTrack_PP }
%
% Colin Glass, Feb 04

handles=guidata(hObject);



[numrows,numcols]=size(handles.MPM);

ImageNamesList = handles.jobvalues.imagenameslist;

dragTailLength = handles.postpro.dragtail;


%starting frame relativ to first frame analysed with polytrack
start = round((handles.postpro.moviefirstimg- handles.jobvalues.firstimage)/handles.jobvalues.increment)+1

%we have to have prior images for dragTail
if start < dragTailLength+2
    start = dragTailLength+2;
end


%last frame relativ to first frame analysed with polytrack
stop = floor((handles.postpro.movielastimg - handles.jobvalues.firstimage)/handles.jobvalues.increment+0.00001)+1 

if stop > handles.jobvalues.lastimage
    stop = handles.jobvalues.lastimage;
end


imagedirectory = handles.jobvalues.imagedirectory;
saveallpath = handles.postpro.saveallpath;

cd(saveallpath)

%initialize the movie
makeQTmovie('start','trackmov.mov')


for movieStep = start:stop
    
    whatcells=[];
    
    %use only the cells chosen by the user, if he chose any
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
     nowImgH=imreadnd2(name,0,handles.jobvalues.intensityMax);

		[rows,cols]= find(handles.MPM);
		
		
		%if the user specified a size for the movie, that's the size we are
		%going to use
		if ~isempty(handles.postpro.figureSize)
		     figure('Position',handles.postpro.figureSize) 
             imshow(nowImgH,[])
        else
             figure, imshow(nowImgH,[])
        end
        
        
        
		hold on
		%one colour per time step (dragTailLength tells you how many
		%timesteps there are
		cmap=jet(dragTailLength+1);

		counter=0;
        
        %loop through the previous pictures (for the tails)
		for i=(2*(movieStep-dragTailLength)):2:(2*movieStep)
            
                counter=counter+1;
                vec=handles.MPM(whatcells,i-3:i);
                [rows,cols]=find(vec==0);
                rows=unique(rows);
                vec(rows,:)=0;
                
                
                
                 for h=1:size(vec,1)
                     
                        if vec(h,1)~=0
							%   plot(vec(h,1),vec(h,2),'.')
							ph=[];
							ph=  plot(vec(h,1:2:3),vec(h,2:2:4));
							set(ph,'Color',cmap(counter,:));
							clear ph
                        end
                      
                end
                
		end
        
        %plot the points that actually belong to the current picture
		plot(vec(:,3),vec(:,4),'r.')
		
		hold off
		
        
        cd(saveallpath)
  
        %add the current figure to the movie
		makeQtmovie('addaxes',gca);
        
		close
end

%finalize the movie
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
% % PROPPERIMG=getinfor(cellprops)
% % 
% % figure,title('coordinates per image')
% % hold on
% % 
% % cellprops

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

