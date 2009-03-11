function [percrange] = boxwhiskerPerRange(tau, k,perc)
% calculate the range of the box-and-whisker plot for a function of
% time constant tau which contains a given percentage of the data
% k indicates if the function is an exponential or a rayleigh


if tau==0
    
    percrange = 0;
    
else
    
    if k==1
        a = 1;
        b = -perc*exp(1);
        c = 1;
        quadsol = ( -b + sqrt( b^2+4 ) )/2;
        factor = log ( quadsol );
        percrange = tau* factor;
    elseif k==2
        if tau>=1
            rvec = [0:0.1:round(3*tau)];
        % if tau<1, use smaller increment
        else
            % exponent
            ep = floor(log10(tau))-1;
            incr = 10^ep;
            rvec = [0:incr:3*tau];
        end
        ivec =  exp(-((tau-rvec).^2)/tau^2) - exp(-((tau+rvec).^2)/tau^2) ;
        rpos = min(find(ivec>=perc));
        percrange = rvec(rpos);
    end

end % of if



end % of function

