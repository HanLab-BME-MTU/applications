function [estimates,renorm,residuals,exFlag,outP,lam,jacobian] = fitFVcurve(param);

options = optimset('Display','off');
start_point = [param.kinesin_probability_matrix(1,2),param.kinesin_probability_matrix(3,1),param.kinesin_probability_matrix(2,1),param.velocity_max_kinesin];
[estimates,renorm,residuals,exFlag,outP,lam,jacobian]= lsqnonlin(@simulateFVcurve, start_point,[],[],options);

    function Difference = simulateFVcurve(params)

        param.kinesin_probability_matrix(1,2) = params(1);
        param.kinesin_probability_matrix(3,1) = params(2);
        param.kinesin_probability_matrix(2,1) = params(3);
        param.velocity_max_kinesin = params(4);

        for F = 1:7;
            param.Force = -F;
            for i = 1:200
                [param,tracks(i)] = positionSimulation(param);
            end
            for i = 1:200
                vel(i) = (tracks(i).points(end,1)-tracks(i).points(1,1))/param.time_max;
            end
            velMean(F) = nanmedian(vel);
            %[velCI(F,:)] = bootci(2000,{bootfun,vel},'alpha',0.05);
        end

        x = 1:7;
        y = -112.5*x+906.25;

        Difference = velMean - y;

    end %of embedded function


param.kinesin_probability_matrix(1,2) = estimates(1);
param.kinesin_probability_matrix(3,1) = estimates(2);
param.kinesin_probability_matrix(2,1) = estimates(3);
param.velocity_max_kinesin = estimates(4);
bootfun = @(x) nanmedian(x);


for p = 1:5;
    param.kinesin_probability_matrix(1,2) = 0.1*p-0.05;

    
    loadParam
for F = 1:7;
    param.Force = -F;
    for i = 1:200
        [param,tracks(i)] = positionSimulation(param);
    end
    for i = 1:200
        vel(i,F) = (tracks(i).points(end,1)-tracks(i).points(1,1))/param.time_max;
    end
    velMean(F) = nanmedian(vel(:,F));
    [velCI(F,:)] = bootci(2000,{bootfun,vel(:,F)},'alpha',0.05);
end

figure, hold on
plot(1:F,velMean,'k')
plot(1:F,velCI(:,1),'k--')
plot(1:F,velCI(:,2),'k--')
x = 1:7;
y = -112.5*x+906.25;
plot(x,y,'r')



end %of function