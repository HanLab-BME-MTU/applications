function [survFunc,pvec,decaytimes] = compSurvivalFunction_general(data1, field1, data2, field2, numBS, typeBS, percvec, id)
% compSurvivalFunction_general compares the survival functions for
% different conditions
%
% SYNOPSIS [survFunc,pvec] = compSurvivalFunction_general(data1, field1,
% data2, field2)
%
% INPUT     data1    = experiment structure with first data set
%           field1   = field in data1 from which the lifetime histogram is
%                      read, e.g. 'survivalFunction_InRegion'
%           data2    = experiment structure with second data set
%           field2   = field in data1 from which the lifetime histogram is
%                      read, E.G. 'survivalFunction_OutRegion'
%           numBS    = number of bootstrap runs (OPTIONAL)
%                      default = 2000
%           typeBS   = type of bootstrap (OPTIONAL)
%                      'movie' (default) = subsample individual movie
%                      combinations
%                      'traj' = subsample individual trajectories (especially
%                      relevant for big difference in inside vs outside)
%           percvec  = vector with percentage points of survival function
%                      for comparison; default is [0.9:-0.1:0.1]
%           id       = set to 1 if you want to compare matched pairs (e.g.
%                      inside/outside of the same movie)
%                   
%                   
%
% OUTPUT    survFunc = survival function 
%                       (first row DATA1, second row Data2)
%           pvec    = vector of p-values;
%                     1: KS-test on difference between data1 and data2
%                     2: 95% threshold for KS-test on bootstrapped
%                     subsamples of data1 vs bootstrapped subsamples of
%                     data2
%                     3. fraction of non-significant boostrap KS-test
%                     p-values 
%           decaytimes  = decay times, i.e. number of frames for which the
%                     respective survival function reaches the percentage
%                     levels specified by percvec; this has the format
%                     row 1: specified percentage levels
%                     row 2: average condition 1
%                     row 3: error condition 1
%                     row 4: average condition 2
%                     row 5: error condition 2
%                     row 6: p-value for t-test between condition 1 and 2
%                     IF ID==1, also added:
%                     row 7: average difference 1-2
%                     row 8: error difference 1-2
%                     row 9: p-value for t-test on difference
%                   
% NOTE: The current version of the function assumes that ALL MOVIES in the
% data structure are acquired at the same framerate, or that the user
% chooses to treat them as being the same, and that the survival functions
% are normalized to the value for the lifetime 1 frame.
%
%
% last modified DATE: 28-Aug-2008 (Dinah)
% last modified DATE: 10-Sep-2008 (Dinah)
% last modified DATE: 02-Oct-2008 (Dinah)


% averaging can only be performed up until the minimum common length of all 
% movies, so we determine the shortest movie length in the structure
n1 = length(data1);
n2 = length(data2);

for i=1:n1
    mlvec1(i) = data1(i).movieLength;
end
minlen1 = min(mlvec1);

for i=1:n2
    mlvec2(i) = data2(i).movieLength;
end
minlen2 = min(mlvec2);

minlen = min(minlen1,minlen2);

idvar = 0;
if nargin>7
    if id==1 
        if n1==n2
            idvar=1;
        else
            disp('idvar can''t be set to 1 since data sets have different sizes');
        end
    end
end
        
        
histmat1 = nan*zeros(length(data1),minlen);
histmat2 = nan*zeros(length(data2), minlen);

% read data for each movie 
for k=1:n1
    
    if isfield(data1,field1)
        survFun = getfield(data1(k), field1);
        histmat1(k,:) = survFun(1:minlen);
    else
        error(['function requires a structure field called ',field1]);
    end          
end

% read data for each movie 
for k=1:n2
    
    if isfield(data2,field2)
        survFun = getfield(data2(k), field2);
        histmat2(k,:) = survFun(1:minlen);
    else
        error(['function requires a structure field called ',field2]);
    end          
end


% average survival functions
survFun1_SUM = nansum(histmat1,1);
survFun2_SUM = nansum(histmat2,1);
survFun1_AVE = survFun1_SUM/n1;
survFun2_AVE = survFun2_SUM/n2;

% define output
survFunc(1,:) = survFun1_AVE;
survFunc(2,:) = survFun2_AVE;


% convert to distributions
df1 = convSurvivalFunction2dist(survFun1_SUM);
df2 = convSurvivalFunction2dist(survFun2_SUM);
nd1 = length(df1);
nd2 = length(df2);

% compare averages
[H,pval_av] = kstest2(df1,df2);




%=================================
%% bootstrap

% number of bootstrap runs
if nargin>4
    nbs = numBS;
else
    nbs = 2000;
end

testtype = 0;
if nargin>5
    if strcmp(typeBS,'traj')
        testtype = 1;
    end
end


% bootstrap first data set
for b=1:nbs
    
    fprintf('bootstrap %03d',round(100*(b/nbs)));
    
    % testtype 0=================================
    % bootstrap from reshuffled movies
    if testtype==0
        % first position vector
        pos1 = randsample(n1,n1,true); 
        % local bootstrap set from normalized matrix
        cmat1 = histmat1(pos1,:);  
        bootstrapSUM1 = nansum(cmat1,1); 
        df1_sub = convSurvivalFunction2dist(bootstrapSUM1);
        % second position vector
        if n2==n1
            pos2 = pos1;
        else
            pos2 = randsample(n2,n2,true); 
        end
        % local bootstrap set from normalized matrix
        cmat2 = histmat2(pos2,:); 
        bootstrapSUM2 = nansum(cmat2,1); 
        df2_sub = convSurvivalFunction2dist(bootstrapSUM2);
        
    else
        
    % testtype 1=================================
    % bootstrap from downsampled trajectories
    
        if nd1>nd2
            % position vector
            pos1 = randsample(nd1,nd2,false);
            df1_sub = df1(pos1);
            df2_sub = df2;
        else
            % position vector
            pos2 = randsample(nd2,nd1,false);
            df1_sub = df1;
            df2_sub = df2(pos2);
        end
        
    end
        
      
    
    %===========================================
    % comp between subsample 1 and total sample 2
    %[H,pval_bs1] = kstest2(df1_sub,df2);
    [H,pval_bs1] = kstest2(df1_sub,df2_sub);
    bootstrap1_pval(b) = pval_bs1;
    
     
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b');
    
end

fprintf('\n');


% sorted bootstrap confidence intervals
sort_bs1 = sort(bootstrap1_pval);

% 95% of all bootstrap comparison p-values are smaller than...
cflevel = sort_bs1(round(0.95*nbs));

% the fraction of bootstrap comparison p-values that are non-significant
% (i.e. larger than 0.05) is...
plevel = length(find(sort_bs1>0.05))/nbs;



% output:
pvec(1) = pval_av; % KS-test p-value of all data1 vs all data 2
pvec(2) = cflevel; % KS-test p-value below which 95% of bootstrap tests are located
pvec(3) = plevel;  % percentage of non-significant boostrap test results

%% level outout

levels = [0.9:-0.1:0.1];
if nargin>6
    if ~isempty(percvec)
        levels = percvec;
    end
end



for k=1:length(levels)
    
    % current percentage increment
    cinc = levels(k);
    %incrVec1 = [];
    %incrVec2 = [];
    % find decay time for this level condition 1
    for a=1:n1
        fpos = find( (histmat1(a,:)/histmat1(a,1))<cinc );
        if ~isempty(fpos)
            incrVec1(a,k) = min( fpos );
        else
            incrVec1(a,k) = nan;
        end
    end
    % find decay time for this level condition 1
    for b=1:n2
        fpos = find( ( histmat2(b,:)/histmat2(b,1) )<cinc );
        if ~isempty(fpos)
            incrVec2(b,k) = min( fpos );
        else
            incrVec2(b,k) = nan;
        end
    end
    
end

for k=1:length(levels)
    
    cincrVec1 = incrVec1(:,k);
    cincrVec2 = incrVec2(:,k);
    
    [h,pval] = ttest2(cincrVec1,cincrVec2);

    decaytimes(1,k) = levels(k);
    decaytimes(2,k) = nanmean(cincrVec1);
    decaytimes(3,k) = nanstd(cincrVec1)/sqrt(n1);
    decaytimes(4,k) = nanmean(cincrVec2);
    decaytimes(5,k) = nanstd(cincrVec2)/sqrt(n2);
    decaytimes(6,k) = pval;
end


%==========================================================================
%% if data are inside vs outside of same movies, compare the individual
%% differences

if idvar==1
    
      
    incrDiff = incrVec1-incrVec2;
    
    for k=1:length(levels)

        cdiff = incrDiff(:,k);
        [h,pval2] = ttest(cdiff); 
        ptdiff(k) = pval2;
        
    end
      
    decaytimes(7,:) = nanmean(incrDiff,1);
    decaytimes(8,:) = nanstd(incrDiff,1);
    decaytimes(9,:) = ptdiff(:);
    
    
    
    figure
    hold on;
    errorbar([1:length(levels)],decaytimes(7,:),decaytimes(8,:),'r.');
    bar([1:length(levels)],decaytimes(7,:),0.3,'r' );

    for n=1:length(levels)
        tlevel(n) = {num2str(levels(n))};
        text( n-0.15,1.1*decaytimes(7,n),['p=',num2str(decaytimes(9,n))] );
    end

    set(gca,'XTick',[1:length(levels)]);
    set(gca,'XTickLabel',tlevel);
    xlabel('percentage level');
    ylabel('decay time pair difference');
    title('single-sided t-test on matched pair differences'); 
    
       
end



%==========================================================================
%% display results

figure; hold on;
cdfplot(df1);
cdfplot(df2);
h = findobj(gca,'type','line');
set(h(1),'linestyle',':','color','r')

xlabel('lifetime (frames)');
ylabel('cumulative fraction');
legend('data1','data2');
axis([0 minlen -0.01 1.01]);
box on


textstr1 = ['ave1 vs ave2: KS-test p=',num2str(pval_av)];

textstr2 = ['bootstrap1 vs bootstrap2 95% p<',num2str(cflevel)];

if length(find(sort_bs1>0.05)) == 0
    textstr3 = ['fraction of non-sign. bootstraps < ',num2str(1/nbs)];
else
    textstr3 = ['fraction of non-sign. bootstraps = ',num2str(plevel)];
end

text(60,0.5,textstr1);
text(60,0.4,textstr2);
text(60,0.3,textstr3);

grid off



figure; hold on;
errorbar([1:length(levels)]-0.15,decaytimes(2,:),decaytimes(3,:),'r.');
errorbar([1:length(levels)]+0.15,decaytimes(4,:),decaytimes(5,:),'b.');
bar([1:length(levels)]-0.15,decaytimes(2,:),0.3,'r' );
bar([1:length(levels)]+0.15,decaytimes(4,:),0.3,'b' );
for n=1:length(levels)
    tlevel(n) = {num2str(levels(n))};
    text( n-0.15,1.1*max(decaytimes(2,n),decaytimes(4,n)),['p=',num2str(decaytimes(6,n))] );
end
set(gca,'XTick',[1:length(levels)]);
set(gca,'XTickLabel',tlevel);
xlabel('percentage level');
ylabel('decay time');
title('double-sided t-test on percentage levels'); 



end % of function




%=========================================================================
%
%                            SUBFUNCTION



function [df] = convSurvivalFunction2dist(sf)
% convSurvivalFunction2dist converts a survival function back into a
% distribution function
%
% SYNOPSIS [df] = convSurvivalFunction2dist(sf)
%
% INPUT     sf    = survival function, presumed at equal spacing 
%
% OUTPUT    df    = distribution function
%
% NOTE: Advantage of this function is that survival functions can be summed
% and cropped until appropriate time points, that that the reconversion
% gets rid of the problem of unequal movie lengths


% convert to histogram
hf = abs(diff(sf));

df = zeros(sum(hf),1);
ct = 0;

for i=1:length(hf)
    numentry = round(hf(i));
    %cv = i + zeros(numentry,1);
    df(ct+1:ct+numentry) = i; 
    ct = ct+numentry;
end % of for

end % of subfunction




    
