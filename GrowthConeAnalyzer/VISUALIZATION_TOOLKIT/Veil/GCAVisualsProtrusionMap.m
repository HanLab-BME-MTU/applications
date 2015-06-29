function [ output_args ] = protrusionMapMine( protSamples,saveDir,makeMovie )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here 
if nargin<3 
    makeMovie =0 ;
end 
%setFigureMoviePlots(200,200,'off')

  sfont = {'FontName','Arial','FontSize',12,'FontName','Arial','FontWeight','bold'};
                lfont = {'FontName','Arial', 'FontSize',14,'FontName','Arial'};
   imgData = protSamples.avgNormal;
                
                
                %subplot(3,2,1:2);
              
                
                % For now just plot yourself but should learn to use ScalarMapDisplay
                % object
                
                Colormap='jet';
                x = colormap;
                x(1,:) = [1,1,1];
                colormap(x)
                
                Colorbar ='on';
                CLim = [-100 100];
                
                % convert imageData to nm per sec
                imgData = imgData.*216; % currently 216 nm per pixel
                imgData = imgData./5;  % 5 sec per frame
               % imgData = imgData(1:130,:); 
   %             imgData = imgData(1:90,:); 
                h= imagesc(imgData,CLim);
            axis xy
                % white out NaN values
                             alphamask =true(size(imgData));
                             alphamask(isnan(imgData))=false;
                             set(h,'AlphaData',alphamask,'AlphaDataMapping','none');
                ylabel('Window Number',lfont{:});
                 xlabel('Time (s)' , lfont{:});
             %  title(name,sfont{:});
               frameNum = get(gca,'XTickLabel');
                 inSec  = str2double(frameNum).*5; % SET THE TIME INTERVAL
                 num2str(inSec);
                 set(gca,'XTickLabel',inSec);
                 set(gca,'FontSize',10,'FontName','Arial');
                 % colorbar
                  
                 
                  if makeMovie~=1
                  saveas(gcf,[saveDir filesep 'WindowsAll.eps'],'psc2');
                  saveas(gcf,[saveDir filesep 'WindowsAll.fig']); 
                  saveas(gcf,[saveDir filesep 'WindowsAll.png']);
                  end % else add scatters 
end

