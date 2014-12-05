function [ h ] = plotNeuriteOutgrowthInTime( neuriteIn ,color,norm,timeInterval,makeMovie,saveDir)
%Overview
%% 
% Notes regarding the Cell: Cell is bad
%% 
if norm == 1 
neuriteIn = neuriteIn-neuriteIn(1); % normalize by the length in the first frame 
end 
%figure('Visible','off'); 

       % scatter((0:length(neuriteIn)-1).*timeInterval,neuriteIn,50,color,'filled'); 
        
        
        scatter((0:length(neuriteIn)-1).*timeInterval,neuriteIn,20,color,'filled'); 
        hold on 
        
       [outgrowthFilt , outlierIdx] = findOutliersFromMedFilt(neuriteIn,6,6);
       
       plot((0:length(outgrowthFilt)-1).*timeInterval,outgrowthFilt,'color','r','Linewidth',2); 
       % plot((0:length(neuriteIn)-1).*timeInterval, neuriteIn,'color',color,'Linewidth',2); 
        xlabel('Time (s)','FontName','Arial','FontSize',20); 
        ylabel({'Net Neurite Outgrowth (um) '} ,'FontName','Arial','FontSize',20); 
        set(gca,'FontName','Arial','FontSize',18); 
        axis([0 (length(neuriteIn)-1)*timeInterval,-10,20]); 
        saveas(gcf,[saveDir filesep 'NetOutgrowth' '.eps'],'psc2'); 
        saveas(gcf,[saveDir filesep 'NetOutgrowth' '.fig']); 
       % close gcf 
        plotVel =0 ; 
       if plotVel == 1
        figure('Visible','off')    
        vel = diff(outgrowthFilt)*1000/5; % should be um/sec
        
        
        plot((0:length(vel)-1).*timeInterval,vel,'color','b','Linewidth',2); 
        ylabel('Velocity (nm/sec)','FontName','Arial','FontSize',20); 
        xlabel('Time (s)', 'FontName','Arial','FontSize',20); 
        saveas(gcf,[saveDir filesep 'Outgrowth_Velocity' '.eps'],'psc2'); 
        saveas(gcf,[saveDir filesep 'Outgrowth_Velociyt' '.fig']); 
       end 
if makeMovie== 1
    
    for iFrame = 1:length(neuriteIn) 
        figure('Visible','off'); 
       
        plot((0:length(neuriteIn)-1).*timeInterval, neuriteIn,'color',color,'Linewidth',2); 
        hold on 
         scatter((iFrame-1)*timeInterval,neuriteIn(iFrame),100,color,'filled'); 
        
        
        xlabel('Time (s)','FontName','Arial','FontSize',12); 
        ylabel({'Net Outgrowth'; 'Difference From First Frame (um)'} ,'FontName','Arial','FontSize',12); 
        axis([0 (length(neuriteIn)-1)*timeInterval, min(neuriteIn(:)),max(neuriteIn(:))]); 
        saveas(gcf,[saveDir filesep 'NetOutgrowth' num2str(iFrame,'%03d') '.png']); 
        close gcf
    end
    
end 
    
    
    
    
    


