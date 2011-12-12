% Spatial sampling functions
%
% Jeffrey L. Werbin 
% Harvard Medical School
%
% Last update: 7/7/2011
%
% This file includes several functions for artifically resampling point
% process data. 
%
% All the functions take the point process as a list of points of 2D
% positions and return a list in the same format.
%

function [list]=PointP_AntibodySimple(StartList,dist, prob)
%PointP_Antibody applies a resampling that mimics antibody staining. This 
%consists of merging points within a certain distance of each other with a
%certain probability. This function is simple because it assumes 100% of
%the protein is bound to an antibody.

%This method is biased. Elements at the begining of the list are more
%likely to be paired then later. Not a problem if point order is not
%spatially biased. 

n = size(StartList,1);

%creates a temporary list and initializes the final list
tList = StartList;
list = [];
i=1;

while(size(tList,1) > 0)
    
    %Chooses point to test, (always the first entry in tList) and removes
    %it from tList
    point = tList(1,:);
    tList = tList(2:end,:);
    
    %finds the distances between StartList(i)  and all other elements
    [temp,d]=knnsearch(point,tList);
    
    %find number of points less than dist
    num = sum((d < dist));
    
    if ( (random('unif',0,1) < prob) && (num > 0) )
        %choose a point to potentially merge with i. All points with distance 
        %of less than dist are equally considered
        [temp,ind]=sort(d);
        partner = ind(random('unid',num));
        
        %if the above test is passed the point i and point point are merged
        %by averaging their position and both point are removed from
        %subsequent searches
        list = [list ; mean([tList(partner,:);point])];
        
        %remove partner from tList
        tList = [tList(1:partner-1,:);tList(partner+1:end,:)];
    else
        list = [list ; point];
    end
end
  
end