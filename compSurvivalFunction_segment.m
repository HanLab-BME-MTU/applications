function [survFunc,tvec] = compSurvivalFunction_segment(data)
% compSurvivalFunction_segment compares the inside and outside survival
% functions in a data structure
%
% SYNOPSIS [survFunc,tvec] = compSurvivalFunction_segment(data)
%
% INPUT     data    = experiment structure, which has to contain the fields
%                       .survivalFunction_InRegion
%                       .survivalFunction_OutRegion
%
% OUTPUT    survFunc = survival function 
%                       (first row IN, second row OUT)
%           tvec    = time vector for survival levels 90-10% (i.e. lifetime
%                       in frames at which 90% survival is reached)
%                       (first row IN, second row OUT)
%
% NOTE: The current version of the function assumes that ALL MOVIES in the
% data structure are acquired at the same framerate, or that the user
% chooses to treat them as being the same, and that the survival functions
% are normalized to the value for the lifetime 1 frame.
%
% To accomodate variable framerates quantitatively, there are two
% possibilities: Either the survival functions can be matched
% pre-averaging, by interpolating the 'slow' survival functions to a higher
% framerate or by binning the 'fast' functions to the slower rate.
% Alternatively, the time vectors can be determined individually for each
% movie in seconds, and then averaged, without averaging the survival
% functions in the first place
%
% last modified DATE: 25-Jun-2008 (Dinah)
% last modified DATE: 23-Jul-2008 (Dinah)


% averaging can only be performed up until the minimum common length of all 
% movies, so we determine the shortest movie length in the structure
for i=1:length(data)
    mlvec(i) = data(i).movieLength;
end
minlen = min(mlvec);


% read survival function for each movie 
for k=1:length(data)
    
    if isfield(data,'survivalFunction_InRegion')
        sfIn = data(k).survivalFunction_InRegion;
        sfOut = data(k).survivalFunction_OutRegion;
    else
        error('function requires existence of a structure field .survivalFunction_InRegion');
    end
    
        
    % normalize survival function
    sfIn_norm = sfIn/sfIn(1);
    sfOut_norm = sfOut/sfOut(1);
    
    sfIn_mat(k,1:minlen) = sfIn_norm(1:minlen);
    sfOut_mat(k,1:1:minlen) = sfOut_norm(1:minlen);
    
end

% average survival function IN
sfIn_av = nanmean(sfIn_mat,1);
sfOut_av = nanmean(sfOut_mat,1);
survFunc(1,:) = sfIn_av;
survFunc(2,:) = sfOut_av;

tauvec = [0.1:0.1:0.9];
tvec = nan*zeros(3,9);
tvec(1,:) = tauvec;

for i=1:length(tauvec)
    pIN = find(sfIn_av<=tauvec(i));
    pOUT = find(sfOut_av<=tauvec(i)); 
    if ~isempty(pIN), tvec(2,i) = min(pIN); end
    if ~isempty(pOUT), tvec(3,i) = min(pOUT); end
end

[H,pval] = kstest2(sfIn_av,sfOut_av);


figure; hold on;
plot(sfIn_av,'r-'); 
plot(sfOut_av,'b-');
plot(tvec(2,:),tauvec,'ro'); 
plot(tvec(3,:),tauvec,'bo');
text(50,0.5,['KS-test pval=',num2str(pval)]);
xlabel('lifetime (frames)');
ylabel('survival function (fraction)');
legend('inside region','outside region');
axis([0 minlen -0.02 1.02]);


% barmat = tvec(2:3,[9:-1:1]);
% for k=1:2
%     barmatDiff(k,2:9) = diff(barmat(k,:));
% end
% barmatDiff(:,1) = barmat(:,1);
% 
% figure; hold on
% barh(barmatDiff,0.4,'stack');
% axis([0 max(barmat(:)) 0.5 2.5]);
% set(gca,'YTick',[1 2]);
% set(gca,'YTickLabel',{'Inside','Outside'});


npt = length(survFunc(1,:));
% calculate base from number of bins
nbin = 12;
bas = exp(log(npt)/(nbin-1));
% calculate number of bins from base
bas = 1.6;
nbin = floor(log(npt)/log(bas))+1;

llogvec = bas.^[0:nbin-1];
ttvec = [1:length(survFunc(1,:))];

ppoints1 = interp1(ttvec,survFunc(1,:),llogvec);
ppoints2 = interp1(ttvec,survFunc(2,:),llogvec);

pdpoints1 = abs(diff(ppoints1));
pdpoints2 = abs(diff(ppoints2));

pdpoints1(length(pdpoints1)+1) = 1-sum(pdpoints1);
pdpoints2(length(pdpoints2)+1) = 1-sum(pdpoints2);

figure
barh([pdpoints1;pdpoints2;(1/nbin)+0*[1:nbin]],0.6,'stack');
axis([0 1 0.4 3.8]);
set(gca,'YTick',[1 2 3]);
set(gca,'YTickLabel',{'Inside','Outside','Colorbar'});
xlabel('fraction in lifetime interval');

for i=1:nbin, names{i}=num2str(round(10*llogvec(i))/10); end
text(((1/nbin)*[0:nbin-1]+0.02),(2.9+0*[1:nbin]),names,'Rotation',90);
text(0.4,3.5,'Lifetime Intervals');

end % of function







    
