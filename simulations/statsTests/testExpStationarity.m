%SCRIPT testExpStationarity generates many trajectories using mtGTPCapLDepK and tests for stationarity

%get experimental data
load wildtype30C;

maxNumSim = length(data);
for i=1:maxNumSim
    totTime(i) = data(i).timePoints(end);
end
totalTime = max(totTime);

mtLengthAve = zeros(totalTime,maxNumSim);
mtLengthSD =  zeros(totalTime,maxNumSim);

for i=1:maxNumSim
    mtLengthAve(data(i).timePoints,i) = data(i).distance(:,1);
    mtLengthSD(data(i).timePoints,i) = data(i).distance(:,2);
end

%save trajectories in file
save('trajectories','mtLengthAve','mtLengthSD');

%calculate the average length and its standard deviation at each time point 
for num = 1:totalTime
    i = find(mtLengthAve(num,:));
    lengthI(num) = length(i);
    mtLengthEnAve(num) = mean(mtLengthAve(num,i)); 
    mtLengthEnSD(num) = std(mtLengthAve(num,i));
end

%calculate the average length and its standard deviation for each sample
for num = 1:maxNumSim
    j = find(mtLengthAve(:,num));
    lengthJ(num) = length(j);
    mtLengthTimeAve(num) = mean(mtLengthAve(j,num));
    mtLengthTimeSD(num) = std(mtLengthAve(j,num));
end

maxShift = 10;  %max tau in autocovarience function

%calculate the autocovarience function using ensemble averaging
autoCov = ones(totalTime,2*maxShift+1)*NaN;
lengthIJ1 = ones(totalTime,2*maxShift+1)*NaN;
for shift = [-maxShift:1:maxShift]
    for t0 = max(1,1-shift):min(totalTime,totalTime-shift)
        t1 = t0 + shift;
        i = find(mtLengthAve(t0,:));
        j = find(mtLengthAve(t1,:));
        ij = intersect(i,j);
        lengthIJ1(t0,maxShift+shift+1) = length(ij);
        autoCov(t0,maxShift+shift+1) = mean((mtLengthAve(t0,ij)-mtLengthEnAve(t0)).*...
            (mtLengthAve(t1,ij)-mtLengthEnAve(t1)));
    end
end

%calculate the autocovarience function using time averaging
ergCov = ones(maxNumSim,2*maxShift+1)*NaN;
lengthIJ2 = ones(maxNumSim,2*maxShift+1)*NaN;
for shift = [-maxShift:1:maxShift]
    for samp = 1:maxNumSim
        tmin = max(1,1-shift);
        tmax = min(totalTime,totalTime-shift);
        i = find(mtLengthAve(tmin:tmax,samp));
        j = find(mtLengthAve(tmin+shift:tmax+shift,samp));
        ij = intersect(i,j)+tmin-1;
        lengthIJ2(samp,maxShift+shift+1) = length(ij);
        ergCov(samp,maxShift+shift+1) = mean((mtLengthAve(ij,samp)-mtLengthTimeAve(samp)).*...
            (mtLengthAve(ij+shift,samp)-mtLengthTimeAve(samp)));
    end
end    
