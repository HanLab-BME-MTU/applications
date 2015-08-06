function [states,kOff,lifetimes] = estimateParameters(receptorInfoLabel,locError)

%ESTIMATEPARAMETERS estimates the binding state and kOff rate for a given set of trajectories
%
%SYNPOSIS [status,kOff,lifetimes] = estimateParameters(receptorInfoLabeled,locError)
%
%INPUT  receptorInfoLabeled: Output of genReceptorInfoLabeled; it is
%                            required to at least have the field
%                            .receptorTraj that contains the position of
%                            each receptor over time. Also, the entries
%                            should be listed in pairs such that the first
%                            one is the red particles from a set, followed
%                            by the green particles, and then repeated for
%                            each video or simulation.
%       locError           : The one-dimensional localization error for all
%                            the particles. More work is needed to change
%                            this for each time point and particle from the
%                            observed data.
%
%OUTPUT states             : A structure array containing the binding
%                            status (1 = monomer, 2 = dimer) for pairs of
%                            red and green particles
%              .status     : This field holds all the statuses for each set
%                            of interesting pairs of red and green
%                            particles. The position in the array is the
%                            sample number it was taken from. The other
%                            matrix dimension is across all the frames
%                            given.
%              .pairs      : Not all possible pairs of red and green
%                            particles are used in the .status field, so
%                            this is a list of which red and which green
%                            particle correspond to which set in the
%                            .status field.
%        kOff              : Returns the estimated kOff rate
%        lifetimes         : A vector with a sorted list of all the
%                            individual lifetimes that were found in the
%                            given data
%
%Code by Paul Blazek June 2015.

%% Parameters

%save space and only keep the receptorTraj field
ril(length(receptorInfoLabel),1) = struct('receptorTraj',[]);
for i = 1:length(receptorInfoLabel)
    ril(i).receptorTraj = receptorInfoLabel(i).receptorTraj;
end
receptorInfoLabeled = ril;

%The number of video/simulation samples given; must be even number of
%entries in the input two get red and green for each sample
numSamples = length(receptorInfoLabeled)/2;
if numSamples-floor(numSamples) ~= 0
    disp('Must have even number of samples');
    return
end

dT = 0.01; %time between frames

%the boxSize is the length of each side of the square in which the
%trajectories are; it is estimated as the maximum position in all of the
%trajectories rounded up to the nearest whole number
boxSize = 0;
for i = 1:numSamples*2
    boxSize = max(boxSize,ceil(max(receptorInfoLabeled(i).receptorTraj(:))));
end

diffConst = 0.1; %diffusion constant (um^2/s)
var1 = 2*locError^2; %variance for state 1 likelihood
var2 = 8*locError^2 + 2*diffConst*dT; %variance for state 2 likelihood
maxState2 = locError*6*sqrt(2); %maximum distance possible for state 2

%The below is not used anymore since it was found to not add any accuracy
%to estimating parameters (in fact it somewhat decreased it)

%this determines the probability that, for the given measured separations,
%the particles are within the cutoff distance
%assocCutoff = 0.05;
%lessThanCutoffProb = zeros(101,2);
%lessThanCutoffProb(:,1) = 0:3*(locError+assocCutoff)/100:3*(locError+assocCutoff); %measured separations
%N = 1500000;
%dist = sqrt((repmat(lessThanCutoffProb(:,1)',N,1)+randn(N,101)*locError*sqrt(2)).^2 + (randn(N,101)*locError*sqrt(2)).^2);
%lessThanCutoffProb(:,2) = mean(dist<assocCutoff);

%% Likelihood Matrices

%these matrices will be used in the estimating parameters step; goodPairs
%are all of the pairs that can possibly bind to each other judged on
%distance
matrices(numSamples,1) = struct('P',[],'piProb',[],'goodPairs',[]);
fprintf('\n');

%this sets up the P and pi matrices for each pair in each sample
for s = 1:numSamples
    progressText(s/numSamples,'Pre-calculations');
    rTraj = receptorInfoLabeled(2*s-1).receptorTraj;
    gTraj = receptorInfoLabeled(2*s).receptorTraj; %
    [numR,~,numFrames] = size(rTraj); %number red receptors
    [numG,~,~] = size(gTraj); %number green receptors
    
    P = zeros(2,2,numFrames,numR,numG,numSamples);
    piProb = zeros(2,numR,numG,numSamples);
    goodPairs = [];
    
    for r = 1:numR
        for g = 1:numG
            %get the separations of the two particles; if they are never close
            %to each other, just skip this pair of particles   
            separations = squeeze(sqrt(sum((rTraj(r,:,:)-gTraj(g,:,:)).^2,2)));
            if sum(separations < maxState2) < 3
                continue
            end
            
            goodPairs = [goodPairs;r g]; %#ok<AGROW>
            
            %set up matrices
            Pi = zeros(2,2,numFrames);
            piProbI = zeros(1,2);
            
            %initial probabilites of being in either state
            piProbI(2) = separations(1)/var1*exp(-separations(1)^2/(2*var1));
            piProbI(1) = 2*separations(1)/boxSize^2;
            
            %calculate P matrix for each time point
            for t = 2:numFrames
                %likelihood for state I
                Pi(2,2,t) = separations(t)/var1*exp(-(separations(t)^2)/2/var1);
                
                %calculate likelihood for state II and scale the matrix P
                if Pi(2,2,t) < 1e-10
                    Pi(2,2,t) = 0;
                    Pi(1,1,t) = 1;
                else
                    for theta = 0:pi/100:2*pi*.9999
                        Pi(1,1,t) = Pi(1,1,t) + exp(-(separations(t)^2 + separations(t-1)^2 - 2*separations(t)*separations(t-1)*sin(theta))/2/var2)*pi/100;
                    end
                    Pi(1,1,t) = Pi(1,1,t)*separations(t)/2/pi/var2;
                    Ps = squeeze(Pi(:,:,t));
                    Pi(:,:,t) = Pi(:,:,t)/geomean(Ps(Ps>0));
                end
            end
            
            P(:,:,:,r,g) = Pi;
            piProb(:,r,g) = piProbI;
        end
    end
    matrices(s).P = P;
    matrices(s).piProb = piProb;
    matrices(s).goodPairs = goodPairs;
    
    clearvars P Pi Ps theta separations goodPairs numR numG rTraj gTraj s
end

%% Estimate Parameters
%iteration of on and off rates to see which set maximizes the
%overall likelihood, using the two-dimensional direct search method
currOn = .5; %current best guess for pOn
currOff = .5; %current best guess for pOff

%upper and lower bounds on pOn and pOff; the scaling is the golden ratio
offUp = 1;
offLo = 0;
onUp = 1;
onLo = 0;
scale = (1 + sqrt(5))/2;

notDone = true;
count = 0;

%run this direct search until the change in parameters is negligible
while notDone
    %reset bounds on pOn and pOff and likelihood
    count = count + 1;
    startOff = currOff;
    startOn = currOn;
    offU = offUp;
    onU = onUp;
    offL = offLo;
    onL = onLo;
    like = 0;
    
    %direct search for best pOff
    progressText(0,sprintf('pOff Run %g',count))
    for i = 1:25
        off1 = offL + (offU-offL)/scale;
        off2 = offU - (offU-offL)/scale;
        like1 = getLikelihood(receptorInfoLabeled,matrices,off1,currOn,maxState2);
        like2 = getLikelihood(receptorInfoLabeled,matrices,off2,currOn,maxState2);
        
        if like1 < like2
            offU = off1;
            like = like2;
        else
            offL = off2;
            like = like1;
        end
        progressText(i/25,sprintf('pOff Run %g',count))
    end
    currOff = (offU+offL)/2; %new best estimate of pOff
    
    %print out temporary result    
    fprintf('currOff = %g - likelihood = %g\n',currOff,like);
    
    %if more work needs to be done, direct search for the best pOn
    if abs(startOff-currOff) >= .01*dT
        progressText(0,sprintf('pOn Run %g',count))
        for i = 1:25
            on1 = onL + (onU-onL)/scale;
            on2 = onU - (onU-onL)/scale;
            like1 = getLikelihood(receptorInfoLabeled,matrices,currOff,on1,maxState2);
            like2 = getLikelihood(receptorInfoLabeled,matrices,currOff,on2,maxState2);
            
            if like1 < like2
                onU = on1;
                like = like2;
            else
                onL = on2;
                like = like1;
            end
            progressText(i/25,sprintf('pOn Run %g\n',count))
        end
        currOn = (onU+onL)/2; %new best estimate of pOn
        
        %print progress
        fprintf('currOn = %g - likelihood = %g\n',currOn,like);
        
        if abs(startOff-currOff) + abs(startOn-currOn) < .02*dT
            notDone = false; %if there wasn't a significant difference, stop
        end
    else
        notDone = false; %if there wasn't a significant difference, stop
    end
end

%display results
fprintf('\nFirst Round:\n');
fprintf('\npOn = %g',currOn);
fprintf('\nkOff = %g\n',currOff/dT);

%% Hidden Markov Model
%get statuses for each good pair and store it in the structure status
states(numSamples,1) = struct('status',[],'pairs',[]);
for s=1:numSamples
    rTraj = receptorInfoLabeled(2*s-1).receptorTraj; %for sample s
    gTraj = receptorInfoLabeled(2*s).receptorTraj; %for sample s
    [~,~,numFrames] = size(rTraj); %get length of video
    goodPairs = matrices(s).goodPairs; %list of relevant pairs of r and g
    
    sampleStatus = zeros(length(goodPairs),numFrames); %set up matrix for statuses
    
    %for each set in the goodPairs matrix, use the hidden Markov model to
    %predict the path of states (1 = monomer, 2 = dimer) for the pairing
    for i = 1:length(goodPairs)
        progressText(i/length(goodPairs),'HMM');
        r = goodPairs(i,1);
        g = goodPairs(i,2);
        sampleStatus(i,:) = hiddenMarkovModel(rTraj(r,:,:),gTraj(g,:,:),locError,currOn,currOff,dT,diffConst,boxSize);
    end
    
    %store this information in the status output structure
    states(s).status = sampleStatus;
    states(s).pairs = goodPairs;
end

[kOff,lifetimes] = getKOff(states,dT); %estimate kOff from the statuses
fprintf('\nRecalculated kOff = %g\n',kOff); %print the estimated off rate





function likelihood = getLikelihood(receptorInfoLabeled,matrices,pOff,pOn,maxState2)

%this calculates the likelihood of seeing everything in receptorInfoLabeled
%given the parameters pOff and pOn
likelihood = 0;
numSamples = length(receptorInfoLabeled)/2;
for s = 1:numSamples
    %get trajectories and matrices for this given video sample
    rTraj = receptorInfoLabeled(2*s-1).receptorTraj;
    gTraj = receptorInfoLabeled(2*s).receptorTraj;
    [~,~,numFrames] = size(rTraj);
    P = matrices(s).P;
    piProb =matrices(s).piProb;
    goodPairs = matrices(s).goodPairs;
    
    for i = 1:length(goodPairs) %loop through the good pairs already found
        r = goodPairs(i,1);
        g = goodPairs(i,2);
        separations = squeeze(sqrt(sum((rTraj(r,:,:)-gTraj(g,:,:)).^2,2)));
        
        %form T matrix, change for each time point to account for
        %separation distances for the on rate, and then
        %incrementally multiply to L, the likelihood matrix
        
        [~,start] = max(separations<maxState2);
        t = max((start-1),2); %the first time point of interest
        if start > 2
            L = squeeze(piProb(:,r,g))'*[1 0; 0 0];
        else
            L = squeeze(piProb(:,r,g))'*eye(2);
        end
        
        pt = zeros(2);
        %Using this while loop proved to speed this up by quite a bit.
        %Basically, a lot of time between pairs is not spent close
        %together, and so this ends up giving a bunch of zeros for the
        %likelihood of being bound. The while loop skips to times when
        %P(1,1) is not zero.
        while t <= numFrames
        %for t = 2:numFrames
            pt(2,1) = P(2,2,t,r,g)*pOff;
            pt(2,2) = P(2,2,t,r,g)*(1-pOff);
            %if (separations(t) < lessThanCutoffProb(end,1))
            %    [~,sepApprox] = min(abs(separations(t)-lessThanCutoffProb(:,1)));
            %    pt(1,2) = P(2,2,t,r,g)*pOn*lessThanCutoffProb(sepApprox,2);
            %else
            %    pt(1,2) = 0;
            %end
            pt(1,2) = P(1,1,t,r,g)*pOn;
            pt(1,1) = P(1,1,t,r,g)-pt(1,2);
            
            L = L*pt;
            
            %this set of statements finds the next t to jump to. If there
            %is no best to jump to, skip to the end
            if t < numFrames
                [maxVal,nextT] = max(separations((t+1):end)<maxState2);
                if maxVal ~= 0
                    if nextT > (t+1)
                        L = L*[1 0; 0 0];
                    end
                    t = t + nextT;
                else 
                    t = numFrames + 1;
                end
            else
                t = numFrames + 1;
            end
            
        end
        
        %get log sum of likelihoods
        L = L*ones(2,1);
        if ~isnan(log(L))
            if isinf(L)
                likelihood = likelihood + log(realmax);
            else
                likelihood = likelihood + log(L);
            end
        end
    end
   
    clearvars rTraj gTraj numFrames P piProb goodPairs r g separations T L
end


function [kOff,allLifetimes] = getKOff(status,dT)

minLifetime = 10; %exclude lifetimes lower than this (the HMM is more 
                  %likely to mess up with these shorter lifetimes
numSamples = length(status);
kOff = zeros(numSamples,1); %vector for the kOff rates
allLifetimes = [];

frames = zeros(numSamples,1);

%% Get Lifetimes

for s = 1:numSamples
    sampleStatus = status(s).status;
    
    %N is the number of tracks, numFrames is the duration of the acquisition
    [N,numFrames] = size(sampleStatus);
    frames(s) = numFrames;
    
    if mean(sampleStatus(:)) == 1 %skip if there are no dimer states
        continue
    end
    
    lifetimes = []; %vector of measured lifetimes
    
    %run over each particle to find the lifetime of all separate binding events
    for n = 1:N
        for t = 1:numFrames %check each time point
            if sampleStatus(n,t) == 2 %the desired status
                life = 1; %start counting the duration
                sampleStatus(n,t) = 1; %change status to prevent double-counting
                
                %a while loop is used to see how long the lifetime is
                keepgoing = true;
                u = t;
                if t == numFrames
                    keepgoing = false;
                end
                while keepgoing
                    u = u + 1;
                    if sampleStatus(n,u) == 2 %if still bound, increase the life
                        life = life + 1;
                        sampleStatus(n,u) = 1;
                        %if it has reached the end, this run won't be used (it
                        %has been found that using these points gives
                        %inaccurate results).
                        if u == numFrames
                            keepgoing = false;
                        end
                    else
                        if t ~= 1
                            lifetimes = [lifetimes; life]; %add lifetime
                        else
                        end
                        keepgoing = false; %stop running because no longer bound
                    end
                end
            end
        end
    end
    cutoffLifetimes = lifetimes(lifetimes>=minLifetime);
    if  isempty(cutoffLifetimes)
        continue
    end
        
    %this sets up the cumulative distribution function
    cdf = zeros(length(cutoffLifetimes),2);
    cdf(:,1) = sort(cutoffLifetimes);
    cdf(:,2) = (1:length(cutoffLifetimes))/length(cutoffLifetimes);
    
    %The equation in this fit is the cumulative distribution of the length of
    %the measured times within a given viewing time frame. t is the tau (or
    %1/kOff) scaled by the dT between frames. v is the total length of the
    %acquisition time. n is the minimum lifetime to include
    myfittype = fittype('(exp(-x/t)*(x+t-v)-exp(-n/t)*(n+t-v))/(exp(-v/t)*t-exp(-n/t)*(n+t-v))','dependent',{'y'},'independent',{'x'},'coefficients',{'t'},'problem',{'v','n'});
    options = fitoptions(myfittype);
    options.Lower = 2;
    options.Upper = (numFrames-1)*5;
    options.StartPoint = 1/dT;
    f = fit(cdf(:,1),cdf(:,2),myfittype,options,'problem',{numFrames,minLifetime});
    
    kOff(s) = 1/(f.t*dT); %store the calculated kOff
    
    allLifetimes = [allLifetimes; lifetimes];
end

allLifetimes = sort(allLifetimes);
cutoffAllLifetimes = allLifetimes(allLifetimes>=minLifetime);

%get the ultimate kOff rate for the whole thing; if they all have the same
%frame length, refit for all the lifetimes, otherwise take the geometric
%mean for all the inidividually calculated kOff's
if numSamples > 1
    if sum(abs(frames(2:end)-frames(1:end-1))) == 0 %if all sample same length
        cdf = zeros(length(cutoffAllLifetimes),2);
        cdf(:,1) = cutoffAllLifetimes;
        cdf(:,2) = (1:length(cdf))/length(cdf);
        f = fit(cdf(:,1),cdf(:,2),myfittype,options,'problem',{numFrames,minLifetime});
        kOff = 1/(f.t*dT);
    else %if the samples have different number of frames
        kOff = geomean(kOff(kOff>0));
    end
end
