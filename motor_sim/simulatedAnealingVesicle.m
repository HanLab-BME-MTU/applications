function [param] = simulatedAnealingVesicle();

timeTotal = 200;
iter = 20;

%MAGNITUDE OF PARAMETER CHANGE
Range12 = 0.1;
Range13 =0.1;
Range31 = 0.1;
Range23 = 0.1;
Range32 = 0.1;
Range22 = 0.1;
RangeMaxVel = 200;
RangeForceStall = 0.5;
RangeStiffness = 0.2;

%LOAD PARAMS
loadParam;

%DATA CURVE TO FIT
F = 1:10;
velData = -112.5*F+906.25;

%CALCULATE FIRST FV CURVE AND ERROR FROM DATA
for f = 1:7;
    param.Force = f;
    for i = 1:iter
        [param,tracks(i)] = positionSimulation(param);
    end
    for i = 1:iter
        vel(i) = (tracks(i).points(end,1)-tracks(i).points(1,1))/param.time_max;
    end
    velMean(f) = nanmedian(vel);
    %[velCI(F,:)] = bootci(2000,{bootfun,vel},'alpha',0.05);
end
output = sum((velMean - velData).^2);

%%
%
for itime = 2:timeTotal
    paramStore = param;
    timePenalty = itime/timeTotal;
    %adjust variables

    changeProb12 = Range12*(rand(1)-0.5);
    if abs(changeProb12/Range12) + timePenalty <= 1
        param.kinesin_probability_matrix(1,2) = param.kinesin_probability_matrix(1,2) +  changeProb12;
    end

    changeProb13 = Range13*(rand(1)-0.5);
    if abs(changeProb13/Range13) + timePenalty <= 1
        param.kinesin_probability_matrix(1,3) = param.kinesin_probability_matrix(1,3) +  changeProb13;
    end

    changeProb31 = Range31*(rand(1)-0.5);
    if abs(changeProb31/Range31) + timePenalty <= 1
        param.kinesin_probability_matrix(3,1) = param.kinesin_probability_matrix(3,1) +  changeProb31;
    end

    changeProb23 = Range23*(rand(1)-0.5);
    if abs(changeProb23/Range23) + timePenalty <= 1
        param.kinesin_probability_matrix(2,3) = param.kinesin_probability_matrix(2,3) +  changeProb23;
    end

    changeProb32 = Range32*(rand(1)-0.5);
    if abs(changeProb32/Range32) + timePenalty <= 1
        param.kinesin_probability_matrix(3,2) = param.kinesin_probability_matrix(3,2) +  changeProb32;
    end

    changeProb22 = Range22*(rand(1)-0.5);
    if abs(changeProb22/Range22) + timePenalty <= 1
        param.kinesin_probability_matrix(2,2) = param.kinesin_probability_matrix(2,2) +  changeProb22;
    end

    changeForceStall = RangeForceStall*(rand(1)-0.5);
    if abs(changeForceStall/RangeForceStall) + timePenalty <= 1
        param.Force_stall_kinesin = param.Force_stall_kinesin + changeForceStall;
    end
    
    changeMaxVel = RangeMaxVel*(rand(1)-0.5);
    if abs(changeMaxVel/RangeMaxVel) + timePenalty <= 1
       param.velocity_max_kinesin = param.velocity_max_kinesin + changeMaxVel;
    end
    
        changeStiff = RangeStiffness*(rand(1)-0.5);
    if abs(changeStiff/RangeStiffness) + timePenalty <= 1
       param.kinesinTetherStiffness = param.kinesinTetherStiffness + changeStiff;
    end


    %SIMULATE FORCE VELOCITY RELATIONSHIP
    tic
    for f = 1:10;
        param.Force = f;
        for i = 1:iter
            [param,tracks(i)] = positionSimulation(param);
        end
        for i = 1:iter
            vel(i) = (tracks(i).points(end,1)-tracks(i).points(1,1))/param.time_max;
        end
        velMean(f) = nanmedian(vel);
        %[velCI(f,:)] = bootci(2000,{bootfun,vel},'alpha',0.05);
    end
    toc

    figure, hold on
    plot(F,velMean,'k')
    plot(F,velData,'r')


    output(itime) = sum((velMean - velData).^2);

    if output(itime) > output(itime-1)
        param = paramStore;
    end %if output is smaller

end %for itime

figure, hold on
plot(F,velMean,'r')
plot(F,velData,'b')
legend('fit','data')
