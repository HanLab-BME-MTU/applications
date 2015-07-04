corr = zeros(25,200);

%set parameters
winSize = 25; %averaging window size
lagT = 1; %lag time for vectors
winLag = winSize+lagT-1; %this is a useful number to have for calculating tail of the shape
numPairs = 1; %number of trajectories to average over

%loop with increasing binding times
for i = 3:25
    actualBoundT = 50; %this is just so each boundT is used twice
    corr(i,:) = troubleshoot(actualBoundT,i,lagT,numPairs); %this gets the corr. coef. from a
    %single run using winSize, lagT, binding time, and the number of runs
    %to average over
    
    %{
    [~,numFrames] = size(corr); %get the number of frames in the simulation
    
    %The following is code cycles through all possible binding time lengths
    %and finds which one fits best
    minError = 5e9;
    for boundT = (winSize+1):numFrames
        
        %initializes shape vector
        shape = ones(winLag+boundT+2,1);
        
        %this creats the general shape of the expected value of the
        %correlation coefficient for a given boundT, winSize, and lagT. The
        %maximum of the shape is set to 1.
        if boundT > winLag
            for j = (1+winLag):-1:(ceil(winLag/2)+2) %left linear
                shape(j) = shape(j+1) - 0.5/floor(winLag/2);
            end
            for j = (boundT+2):(boundT+1+ceil(winLag/2)) % right linear
                shape(j) = shape(j-1) - 0.5/ceil(winLag/2);
            end
            for j = (ceil(winLag/2)+1):-1:2 %left exponential
                shape(j) = shape(j+1)*5/7;
            end
            for j = (boundT+2+ceil(winLag/2)):(length(shape)-1) %right exponential
                shape(j) = shape(j-1)*5/7;
            end
            shape(1) = 0; %starting
            shape(end) = 0; %ending
        elseif boundT >= ceil(winLag/2)
            
        else
            
        end
        
        
        for coeff = 0.5:0.01:1 %for different heights of the shape to maximize fit
            newShape = shape*coeff; %scale the shape
            for offset = (winLag+2):(numFrames-length(newShape)+1) %offset the shape across length
                %this is the squared error
                error = mean((newShape-corr(k,offset:(offset+length(newShape)-1))').^2);
                if error < minError %if this is the best fit so far, save values
                    bestCoeff = coeff;
                    bestOffset = offset;
                    bestBoundT = boundT;
                    minError = error;
                    bestShape = newShape;
                end
            end
        end
    end
    
    startT = bestOffset + ceil(winLag/2)+1; %given the different offsets, this is the start time
    
    %print the results
    fprintf('\nFound boundT of %d for actual %d',bestBoundT,actualBoundT)
    fprintf('\nFound startT of %d for actual %d',startT,round((numFrames-actualBoundT)/2))
    fprintf('\nMean error: %g\n',minError)
    
    %for future use: this subtracts off the shape from the correlation
    %graph so that subsequent passes can be made over the data
    corr(k,bestOffset:(bestOffset+length(bestShape)-1)) = corr(k,bestOffset:(bestOffset+length(bestShape)-1)) - bestShape';
    corr(k,:) = corr(k,:).*(corr(k,:)>0);
    %}
end

clear minError winSize lagT winLag shape boundT newShape coeff offset minError numFrames
