function corr = distanceCorrelation(sample1,sample2,time,winSize)

sample1 = squeeze(sample1);
sample2 = squeeze(sample2);

%checks that x and y and the same size
if size(sample1) ~= size(sample2)
    disp('Inputs must be same size')
    return
end

[length, dim] = size(sample1);
if length < dim
    sample1 = sample1';
    sample2 = sample2';
end

sample1 = sample1(time:(time+winSize-1),:);
sample2 = sample2(time:(time+winSize-1),:);

corr = covariance(sample1,sample2)/sqrt(covariance(sample1,sample1)*covariance(sample2,sample2));



function cov = covariance (sample1,sample2)
[length, ~] = size(sample1);

a = zeros(length);
b = a;

for i=1:length
    a(i,:) = sqrt((sample1(i,1)-sample1(:,1)).^2 + (sample1(i,2)-sample1(:,2)).^2);
    b(i,:) = sqrt((sample2(i,1)-sample2(:,1)).^2 + (sample2(i,2)-sample2(:,2)).^2);
end

A = a - (repmat(mean(a,1),length,1) + repmat(mean(a,2),1,length) - mean(a(:)));
B = b - (repmat(mean(b,1),length,1) + repmat(mean(b,2),1,length) - mean(b(:)));

cov = sqrt(sum(sum(A.*B)))/length;