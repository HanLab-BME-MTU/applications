function [p]=ttestBootstrap(pop1,pop2)
%close all

% bootstrapped t-test

nReps=1000;

n1=length(pop1);
n2=length(pop2);

u1=zeros(nReps,1);       u2=zeros(nReps,1);
std1=zeros(nReps,1);     std2=zeros(nReps,1);
skew1=zeros(nReps,1);    skew2=zeros(nReps,1);
for i=1:nReps
    % sample both populations with replacement
    s1=pop1(ceil((n1-1).*rand(n1,1)));
    s2=pop2(ceil((n2-1).*rand(n2,1)));
    
%    s1=randsample(pop1,n1,true);
 %   s2=randsample(pop2,n2,true);

    % get mean, std, and skewness of both populations
    u1(i)=mean(s1);
    u2(i)=mean(s2);
    std1(i)=std(s1);
    std2(i)=std(s2);
    skew1(i)=skewness(s1);
    skew2(i)=skewness(s2);

end

% m1=mean(u1);
% m2=mean(u2);
% z1=std(u1);
% z2=std(u2);
% 
% pop3=m1+z1.*randn(100,1); figure; hist(pop3,20);
% pop4=m2+z2.*randn(100,1); figure; hist(pop4,20);
% [h,p1]=ttest2(pop3,pop4,0.05,'both');

for iParam=1:3
    switch iParam
        case 1
            test='mean';
            stat1=u1;
            stat2=u2;
        case 2
            test='std';
            stat1=std1;
            stat2=std2;
        case 3
            test='skewness';
            stat1=skew1;
            stat2=skew2;
    end

    h1 = lillietest(stat1);
    h2 = lillietest(stat2);

    if h1==0 && h2==0
        [h,p1]=ttest2(stat1,stat2,0.05,'both','unequal');
    else
        p1=nan;
        disp(['One or both ' test ' populations failed the Lilliefors test for normality']);
    end

    switch iParam
        case 1
            p.pMeanBoot=p1;
        case 2
            p.pStdBoot=p1;
        case 3
            p.pSkewBoot=p1;
    end

end



% PERMUTATION TEST

% get union of the two populations
bigPop=[pop1;pop2];

% get absolute value of the difference between the actual population means
deltaPop=abs(mean(pop1)-mean(pop2));

delta=zeros(nReps,1);
for i=1:nReps
    s1=randsample(bigPop,n1,'true');
    s2=randsample(bigPop,n2,'true');

    delta(i)=mean(s1)-mean(s2);

end

p.pMeanPerm=1-normcdf(deltaPop,mean(delta),std(delta));





% DISTRIBUTION TEST

p1=zeros(nReps,1);       
p2=zeros(nReps,1);
for i=1:nReps
    % sample both populations with replacement, but sample the first twice
    s1a=randsample(pop1,n1,true);
    s1b=randsample(pop1,n1,true);
    
    s2a=randsample(pop2,n2,true);
    
    [h p1(i)]=kstest2(s1a-mean(s1a),s1b-mean(s1b)); % calibration
    [h p2(i)]=kstest2(s1a-mean(s1a),s2a-mean(s2a));
end
% figure; hist(p1,20)
% figure; hist(p2,20)
thresh=prctile(p1,5);
p.confDistribDiff=100*(sum(p2<thresh)/length(p2));


subplot(2,6,1); hist(pop1,20); xlabel('pop1')
subplot(2,6,7); hist(pop2,20); xlabel('pop2')
subplot(2,6,2); hist(u1,20); xlabel('pop1 mean')
subplot(2,6,8); hist(u2,20); xlabel('pop2 mean')
subplot(2,6,3); hist(std1,20); xlabel('pop1 std')
subplot(2,6,9); hist(std2,20); xlabel('pop2 std')
subplot(2,6,4); hist(skew1,20); xlabel('pop1 skew')
subplot(2,6,10); hist(skew2,20); xlabel('pop2 skew')
subplot(2,6,5); hist(delta,20); xlabel('perm delta')
%subplot(2,6,11); hist(s2,20); xlabel('pop2 perm')
%subplot(2,6,6); hist(p1,20); xlabel('dist calib')
%subplot(2,6,12); hist(p2,20); xlabel('dist delta')

