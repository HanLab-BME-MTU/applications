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

function [list]=PointP_Antibody(StartList, eff, dist, prob)
%PointP_Antibody applies a resampling that mimics antibody staining. This 
%consists of picking a subset of the points that an antibody binds and then
%determines if it binds a second point within a certain distance with a
%certain probability

n = size(StartList,1);

%tList is the initial list of molecules with a bound antibody
%pList is the list of remaining molecules and list is the output
tList = PointP_SubSample(StartList,eff);
pList = setdiff(StartList,tList,'rows');
list = [];

while(size(tList,1) > 0)
    
    %Chooses point to test, (always the first entry in tList) and removes
    %it from tList
    point = tList(1,:);
    tList = tList(2:end,:);
    
    %finds the distances between StartList(i)  and all other elements
    [temp,d]=knnsearch(point,pList);
    
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
        list = [list ; mean([pList(partner,:);point])];
        
        %remove partner from tList
        pList = [pList(1:partner-1,:);pList(partner+1:end,:)];
    else
        list = [list ; point];
    end
end
    
end