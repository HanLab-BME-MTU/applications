function speed(hObject)

handles=guidata(hObject);

bins=[0:1:60];

[numrows,numcols]=size(handles.MPM);

dist=zeros(numrows,1);

saveallpath=get(handles.GUI_fm_saveallpath_ed,'String');

for i=1:2:numcols-3
        
        vec=handles.MPM(:,i:i+3);
        [rows,cols]=find(vec==0);
        rows=unique(rows);
        vec(rows,:)=0;
        
        numbercells((i+1)/2,1)=(numrows-length(rows));
        dist(:,(i+1)/2)=sqrt((vec(:,1)-vec(:,3)).^2 + (vec(:,2)-vec(:,4)).^2);
        avardist((i+1)/2,1)=sum(dist(1:numrows,(i+1)/2))/(numrows-length(rows));
        
        realdist=dist(find(dist(:,(i+1)/2)),(i+1)/2);
        speedvar((i+1)/2)=var(realdist);
         
        velHist(:,(i+1)/2)=hist(realdist,bins)';
        velHist(:,(i+1)/2)=velHist(:,(i+1)/2)./(numrows-length(rows));
    end
       
    
        h_fig=figure
        surf(velHist)
        hgsave(h_fig,[saveallpath filesep '3dSpeedHistBar.fig']);
		print(h_fig, [saveallpath filesep '3dSpeedHistBar.eps'],'-depsc2','-tiff');
		print(h_fig, [saveallpath filesep '3dSpeedHistBar.tif'],'-dtiff');
		
       
        h_fig=figure
        bar3(bins,velHist)
		hgsave(h_fig,[saveallpath filesep '3dSpeedHistBar.fig']);
		print(h_fig, [saveallpath filesep '3dSpeedHistBar.eps'],'-depsc2','-tiff');
		print(h_fig, [saveallpath filesep '3dSpeedHistBar.tif'],'-dtiff');
		
        
        
        
		h_fig=figure,plot(numbercells),title('number of cells')
		xlabel('Frames')
		ylabel('# cells')
		ymax=max(numbercells);
		axis([0 numcols/2 0 ymax])
		
		
		hgsave(h_fig,[saveallpath filesep 'numbercells.fig']);
		print(h_fig, [saveallpath filesep 'numbercells.eps'],'-depsc2','-tiff');
		print(h_fig, [saveallpath filesep 'numbercells.tif'],'-dtiff');
		
		
		h_fig=figure
		subplot(2,1,1); plot(avardist),title('avarage velocity of cells')
		xlabel('Frames')
		ylabel('distance per frame (in pixel)')
		ymax=max(avardist);
		axis([0 numcols/2 0 ymax])
		
	
		subplot(2,1,2); plot(speedvar),title('variance of avarage velocity of cells')
		xlabel('Frames')
		ylabel('variance')
		ymax=max(speedvar);
		axis([0 numcols/2 0 ymax])
		
        
		hgsave(h_fig,[saveallpath filesep 'velocity.fig']);
		print(h_fig, [saveallpath filesep 'velocity.eps'],'-depsc2','-tiff');
		print(h_fig, [saveallpath filesep 'velocity.tif'],'-dtiff');
		
    
   
    
   
cd(saveallpath)
save('distances', 'dist')
save('speedvariances','speedvar')
save('avardist','avardist')
    
   guidata(hObject, handles);