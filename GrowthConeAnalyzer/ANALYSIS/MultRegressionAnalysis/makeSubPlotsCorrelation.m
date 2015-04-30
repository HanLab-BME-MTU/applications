function [ output_args ] = makeSubPlotsCorrelation(output)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% get all p-Values 
pValues = output(:).P; 

% filter 
corr = output(pValues<0.05); 

potentialCorr = output(pValues>0.05 & pValues <0.01);

noCorr = output(pValues>0.1); 

mFont= {'FontSize',20,'FontType','Arial'}; 

 nParams = length(corr); 
nRows = ceil(nParams/3) ;
add = 0 ;
fsfigure(1);
for i = 1:nRows
    
subplot(nRows,3,i+add); 
scatter(output(i).descriptor,output(i).response,'filled','k'); 
xlabel(output(i).name,mFont{:}); 
ylabel('Mean Velocity Per Growth Segment (um/min)',mFont{:}); 
end 
save(gcf,'corr.fig'); 
end 



