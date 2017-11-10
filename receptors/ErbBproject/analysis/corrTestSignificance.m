function [ correlation,confidence ] = corrTestSignificance( data,varargin)
%corrTestSignificance takes in data used for correlation anaylsis and
%determines signficance level as compared to random shuffling of pairs
%   Detailed explanation goes here

  ip=inputParser;
ip.CaseSensitive=false;
ip.StructExpand=true;

ip.addRequired('data',@isnumeric);
ip.addOptional('alpha',0.01,@isnumeric);
ip.addOptional('Robust',false,@islogical);

ip.parse(data,varargin{:});

% supress warnings
w = warning ('off','all');

alpha = ip.Results.alpha;
Robust= ip.Results.Robust;

n = size(data,1);
h = size(data,2);

if Robust
corrTemp = @corrRobust;
 else
   corrTemp = @corr;
end

dist = zeros([h,h,500]);

for i=1:500
    for j=1:h
    ind=randperm(n);
    temp = data;
    temp(1:n,j)=data(ind,j);
    
    c = corrTemp(temp);
    dist(j,:,i)= c(j,:);
    end
end

%find signficance by ordering the values for the perm corr and find the
%alpha/2 extremes these are the confidence intervals
sig = floor((alpha/2)*500);

confidence = zeros([h,h,2]); %upper bound is 1 lower bound is 2
for i=1:h
    for j=1:h
        s= sort(dist(i,j,:));
        confidence(i,j,1) = s(1+sig);
        confidence(i,j,2) = s(end-sig);
    end
end
    
correlation = corrTemp(data);

w = warning ('on','all');

end


function c= corrRobust(data)

s = size(data);
n=s(2);
c = zeros(n);

for i=1:n
	for j=i:n
		[b,stats]=robustfit(data(:,i),data(:,j));
                c(i,j)= b(2);
                c(j,i)= b(2);
        end
end


end
