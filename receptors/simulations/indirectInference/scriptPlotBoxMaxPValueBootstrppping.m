%script to boxplot maxPvalue in funcion of rD,aP and lR

%%
currDir ='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170112/bootstrapping/results';
% the name until target
% title of the figure will consider the infos that are filled here for the
% name of the target
 rDtarget = {'rD10'};%,,'rD40','rD60','rD80','rD100','rD120','rD140','rD160'};
 aPtarget = {'aP0p5'};%,'aP0p4','aP0p5','aP0p6','aP0p7','aP0p8'};
 lRtarget ={'lR0p3'};%,'lR0p3','lR0p4','lR0p5'};


%% plot

% Call cell with max p-value and respective coordinates (only significant values)
temp= load([currDir,filesep,rDtarget{1},aPtarget{1},lRtarget{1 },filesep,'maxPvalueCoordSig.mat']);
maxPvalueCoor=temp.maxPvalueCoordSig;

%separate p-value and coordinates into two vectors
MaxPvalues=vertcat(maxPvalueCoor{:,1});
coordMaxPvalue=vertcat(maxPvalueCoor{:,2});

%counts for coordinates that appears more than 10 times
[uniqueRows,~,ind] = unique(coordMaxPvalue,'rows');
counts = histc(ind,1:max(ind));

%replace by NaN the rows with less than 10 repetitions
condition=find(counts<10);

if condition~=0
for indexR=1:length(condition)
ind(ind==condition(indexR))=NaN;   
end
end


%infoBoxPlot
infoBoxPlot=[uniqueRows,counts];
%boxPlot
boxplot(MaxPvalues,ind)
ylabel('p-value','FontSize',22); 
xlabel('[rD,aP,lR]','FontSize',22)
set(gca,'FontSize',22)
% title(['target:',rDtarget{1},filesep,aPtarget{1},filesep,lRtarget{1}],'FontSize',12);
figH = gcf; 

      
 %% Save figure

        outFile = [currDir,filesep,rDtarget{1},aPtarget{1},lRtarget{1 },filesep,'maxPvalueBoxPlot'];
        saveas(figH,outFile,'png');
        saveas(figH,outFile,'fig');  
        
  % save file with the information of equivalent coordinates for the box plot
  resultsDirBoot=[currDir,filesep,rDtarget{1},aPtarget{1},lRtarget{1}];
   save([resultsDirBoot,filesep,'infoBoxPlot'],'infoBoxPlot','-v7.3');  
 