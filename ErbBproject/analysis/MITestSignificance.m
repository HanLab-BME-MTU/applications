function [ MI,confidence ] = MITestSignificance( data,alpha )
%MITestSignificance takes in data used for correlation anaylsis and
%determines signficance level as compared to random shuffling of pairs
%   Detailed explanation goes here

n = size(data,1);
h = size(data,2);

dist = zeros([h,h,500]);

for i=1:500
    for j=1:h-1
        ind=randperm(n);
        A = data(ind,j);
        for k = j+1:h
            dist(j,k,i)= mutualInfo(A,data(:,k));
            dist(k,j,i)= dist(j,k,i);
        end
    end
end

%find signficance by ordering the values for the perm corr and find the
%alpha/2 extremes these are the confidence intervals
sig = floor((alpha/2)*500);

confidence = zeros([h,h,2]); %upper bound is 1 lower bound is 2
MI = zeros([h,h]);
for i=1:h
    for j=1:h
        %calculate MI
        MI(i,j)= mutualInfo(data(:,1),data(:,j));
        
        %calculate confidence
        s= sort(dist(i,j,:));
        confidence(i,j,1) = s(1+sig);
        confidence(i,j,2) = s(end-sig);
    end
end
    MI = corr(data);
end

function [MI] =mutualInfo(A,B)
%calculates mutual Information between A and B
%the joint distribution is calculated as a 50x50 histogram with bins
%ranging from 0 and the value of the 90th precentile

%find ranges for A and B
s = sort(A);
n = numel(A);
maxA = s(floor(n*.9));
rangeA = [0:maxA/50:maxA];


s = sort(B);
n = numel(B);
maxB = s(floor(n*.9));
rangeB = [0:maxB/50:maxB];

[joint,bin]= hist3([A,B],'Edges',{rangeA,rangeB});

MI = sum(sum(joint.*log2(joint)));

end
