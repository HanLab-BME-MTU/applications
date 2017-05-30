% script to save all the info from the simulations

currDir='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170501/probe';

% cellInfoAllProbes=cell(100,5);
% count to know which row will be filled
% temp = load(fullfile(currDir,filesep,'cellInfoAllIntStatic.mat'));
% cellInfoAllIntStatic=temp.cellInfoAllIntStatic;
%    indexCount=1+size(cellInfoAllSim,1);
%  indexCount=5;
%fill the columns with the info

%first column is the path
% cellInfoAllSim{1,1}='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170327/varyDissRate/targetISruns';
%  cellInfoAllSim{indexCount,1}='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170424/targetISruns';
% % 
% % % %second column is the receptor density
%  cellInfoAllSim{indexCount,2}={'rD14'};
%  
%  cellInfoAllSim{indexCount,6}=14;
% 
% % %
% % % %third column is the association probability
% % %
%  cellInfoAllSim{indexCount,3}={'aP0p3','aP0p7'};%'aP0p4','aP0p6'
%  
%  cellInfoAllSim{indexCount,7}=[0.3,0.7];%0.4,0.6
% % %
% % % %fourth column is the value of dissociation rate
%  cellInfoAllSim{indexCount,4}={'dR1'};
%  
%  cellInfoAllSim{indexCount,8}=1;
% % %
% % % %fifth column is the value of labeled fraction
%  cellInfoAllSim{indexCount,5}={'lR0p2','lR0p4'};
% 
% 
% % %third column is the association probability
% % 
%  
% % 
% % %fourth column is the value of dissociation rate
%  
% % 
% % %fifth column is the value of labeled fraction
%  cellInfoAllSim{indexCount,9}=[0.2,0.4];
% % 
% % 
% % % outPut numbers
%  cellInfoAllSim{indexCount,10}=1:10;

cellInfoAllIntermidiateStatatistics=cellInfoAllIntermidiateStatatistics(:,[1 2 4 3 5 6 7 8 9 10]);
%save updated cell array

save([currDir,filesep,'cellInfoAllIntermidiateStatatistics'],'cellInfoAllIntermidiateStatatistics','-v7.3');