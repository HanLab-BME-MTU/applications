function [omega,omegaHigh,omegaLow]=brInvShrinkSpline(velData,velDiscret,omegaDiscret);
% brInvShrinkSpline return omega for MT depoly. according to peskin & oster 1995



maxIndice=find(velDiscret==max(velDiscret));



% 2 main case monotonic or not
%if monotonic, velmax @ 1
if maxIndice==1
    omegaLow=NaN;
    omegaHigh=NaN;
    
    % interpolation thanks spline
    polynome=spline(velDiscret,omegaDiscret);
    omega=ppval(polynome,velData);

% if non monotonic
elseif maxIndice>1
    omegaHigh=ones(size(velData))*NaN;
    omegaLow =ones(size(velData))*NaN;
    omega=NaN;
    indece2Omega=find(velData>velDiscret(1));
    
    if length(indece2Omega)>1
        polynomeLow=spline(velDiscret(1:maxIndice),omegaDiscret(1:maxIndice));
        omegaLow(indece2Omega)=ppval(polynomeLow,velData(indece2Omega));
    end
    %DEBUG
    if length(velDiscret(maxIndice:end))<3
        omegHigh=NaN;
        omegaLow=NaN;
        polynome=spline(velDiscret,omegaDiscret);
        omega=ppval(polynome,velData);
    else
        polynomeHigh=spline(velDiscret(maxIndice:end),omegaDiscret(maxIndice:end));  
        omegaHigh=ppval(polynomeHigh,velData);
    end
end

    
    
    
    




