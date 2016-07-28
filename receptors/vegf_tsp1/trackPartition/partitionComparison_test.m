[CML1,path1] = uigetfile('*.mat','Select a CML or ML');
% [path2,CML2] = uigetfile('*.mat','Select a CML or ML');
partitionComparison('noVEGF',[path1,CML1]);
