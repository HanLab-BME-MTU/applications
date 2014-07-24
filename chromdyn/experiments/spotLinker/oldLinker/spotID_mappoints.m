function [out]=spotID_mappoints(spots,nspots,samp,sxyz,deltamp,weight)
%spotID_mappoints maps corresponding spots by finding the permutation with least difference in distance and intensities 
%
% SYNOPSIS  [idmap.spot.linkdown, idmap.spot.linkup]=spotID_mappoints(spots,nspots,sigma,deltamp)
%
% INPUT spots : list of coordinates and amplitudes of spots of two consecutive frames handed over by spotID
%       nspots : 1*2 vector with number of spots in respective frames
%       sigma.xyz : std of difference in coordinates between spots of two consecutive frames
%       sigma.amp : std of difference in intensities between spots of two consecutive frames
%       deltamp   : total amplitude decrease from t1 to t2
%       weight    : relative weight of distance and intensiy differences (weight distance/weight intensity; w=[0...1...inf])
%
% OUTPUT out.linkdown : mapping of points of frame 1 to frame 2
%        out.linkup   : mapping of points of frame 2 to frame 1
%
% DEPENDS ON 
%
% CALLED BY spotID
%
% SUBFUNCTIONS generatePerm_LuT, spotweight
%   
% c: 23/09/02	Jonas


%-------------initialize

sigma.amp=samp;
sigma.xyz=sxyz;

switch nspots(1)<nspots(2) %if number of spots increases, swap indices 
case 0
    swap=0;
case 1
    swap=1;
    temp.amp=spots(1).amp;
    temp.xyz=spots(1).xyz;
    spots(1).amp=spots(2).amp;
    spots(1).xyz=spots(2).xyz;
    spots(2).amp=temp.amp;
    spots(2).xyz=temp.xyz;
    nspots=nspots'*[0,1;1,0];
end

deltan=nspots(1)-nspots(2);

if deltan<0
    error('wrong swap or similar problem')
end

%generate look-up-table of permutations
perm_LuT=generatePerm_LuT(nspots(2),deltan);

%generate distance matrix
dxyzmatrix=distMat2(spots(1).xyz,spots(2).xyz);

%-------------calculate distances between points

for i=1:size(perm_LuT,1) %for each permutation
    idx2=perm_LuT(i,:); %idx2(2)=index of spot in second frame onto which spot 2 in first frame is mapped
    %% not used %% idxDistMat=sub2ind([nspots(1) nspots(2)],1:nspots(1),perm_Lut(i,:))'; %corresponding indices in distance matrix
    wij=spotweight(spots,dxyzmatrix,sigma,deltamp,idx2,weight);
    
    distance(i)=wij/length(idx2);%normalize by number of links to have it comparable
        %old version (dT): mdist(i)=sum((distMat(matInd)./wij').^2);
end;

[mindist mindistidx]=min(distance);

%------------write output
out.swap=swap;
out.pL=perm_LuT(mindistidx,:);
out.dist=mindist;







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function perm_LuT=generatePerm_LuT(n2,dn)
%generates permutation look-up-table 


%initialize strings
str1='';
str2='PermLut=[perms([1:n2]); PermLut];';
str3='';
PermLut=[];

%for every disappearing spot X->Y: add one loop and one additional variable
% with dn=1: for i1=1:n2;PermLut=[perms([1:n2,i1]);PermLut];end
for i=1:dn
    str1=[str1,'for i',num2str(i),'=i',num2str(i-1),':n2;'];
    str2=[str2(1:end-13),',i',num2str(i),']); PermLut];'];
    str3=[str3,'end;'];
end

%evaluate function to generate Permutation Look-up Table
i0=1;
eval([str1,str2,str3]);
perm_LuT= unique(PermLut,'rows');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function w=spotweight(spots,dxyz,sigma,deltamp,idx2,weight);
%SPWEIGHT compute vector of weights for current configuration

%c: 20/08/02 dT


sumAmp=0;
sumXyz=0;


for i=1:length(idx2)
    nspot=length(find(idx2(i)==idx2)); %how many spots point to spot idx2(i)
    
    %difference in amplitudes= 'amp of spot @t2'-'decrease of intensity'-'amps of spots mapped onto spot @t2'
    damp=abs((spots(2).amp(idx2(i))-deltamp)-sum([spots(1).amp([find(idx2(i)==idx2)])]));
    %weight of permeation/configuration i=exp(-(delta amp/sigma amp)^2)*exp(-(delta xyz/sigma xyz)^2), weighted with weight
    %best configuration has the lowest distance w
    sumAmp=sumAmp+damp^2;
    sumXyz=sumXyz+dxyz(i,idx2(i))^2;
    
end;

w=1-(exp(-(sumAmp*(1-weight))/(sigma.amp)^2))*(exp(-(sumXyz*weight)/(sigma.xyz)^2));