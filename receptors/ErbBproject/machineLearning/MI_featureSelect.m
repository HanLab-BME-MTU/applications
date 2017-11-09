function [ MI,rand,rand_std ] = MI_featureSelect(data,type )
%MI_featureSelect takes in data used for correlation anaylsis and
%determines signficance level as compared to random shuffling of pairs
%
% type must contain only integers

n = size(data,1);
h = size(data,2);

classes = unique(type);
c = numel(classes);

dist = zeros([h,500]);

%for i=1:500
%    for j=1:h
%        ind=randperm(n);
%        A = data(ind,j);
%        dist(j,i)= mutualInfo(A,type);
%    end
%end

rand = mean(dist,2);
rand_std = sqrt(var(dist,[],2));


MI = zeros([h]);
for i=1:h
        %calculate MI
        MI(i)= mutualInfo(data(:,i),type);    
end

 %   MI = corr(data);
end

function [MI] =mutualInfo(A,type)
%calculates mutual Information between A and B
%the joint distribution is calculated as a 50x num Classes histogram with bins
%ranging from 0 and the value of the 90th precentile

%find ranges for A and B
s = sort(A);
n = numel(A);
maxA = s(floor(n*.9));
rangeA = [0:maxA/50:maxA];


rangeB = unique(type);
%rangeB = [0,s+0.1];

[joint,bin]= hist3([A,type],'Ctrs',{rangeA,rangeB});
joint=joint/sum(joint(:));

MI = joint.*log2(joint./(sum(joint,2)*sum(joint,1)));
MI(isnan(MI))=0;
MI = sum(sum(MI));

end
