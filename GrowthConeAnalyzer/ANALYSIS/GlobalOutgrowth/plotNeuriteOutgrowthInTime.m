function [ output_args ] = plotNeuriteOutgrowthInTime( neuriteIn ,color,norm,timeInterval,makeMovie,saveDir)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if norm == 1 
neuriteIn = neuriteIn-neuriteIn(1); % normalize by the length in the first frame 
end 
%figure('Visible','off'); 

       % scatter((0:length(neuriteIn)-1).*timeInterval,neuriteIn,50,color,'filled'); 
        hold on
        scatter((0:length(neuriteIn)-1).*timeInterval,neuriteIn,50,color,'filled'); 
        plot((0:length(neuriteIn)-1).*timeInterval, neuriteIn,'color',color,'Linewidth',2); 
        xlabel('Time (s)','FontName','Arial','FontSize',12); 
        ylabel({'Net Neurite Outgrowth'; 'Normalized To First Frame'} ,'FontName','Arial','FontSize',12); 
        axis([0 (length(neuriteIn)-1)*timeInterval, min(neuriteIn(:)),max(neuriteIn(:))]); 
        saveas(gcf,[saveDir filesep 'NetOutgrowth' '.eps'],'psc2'); 
        saveas(gcf,[saveDir filesep 'NetOutgrowth' '.fig']); 
        close gcf 
        
        
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
    
    
    
    
    


