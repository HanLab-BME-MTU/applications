function koff = koff_unified(f,fThresh)
% Calculates koff as a function of force for a5b1 integrins (kong et al.
% JCB). Force in N, koff in s-1.

    if nargin<2
        fThresh = 0.5e-13; % Default for arp2/3 inhibition
    end

    f = abs(f);

    a = 0.*f;
    b = 0.*f;    
    
    a(f <= fThresh) = 0.1905;
    b(f <= fThresh) = 2.333e+010;
    
    a(f > fThresh) = 0.1905;
    b(f > fThresh) = 7e+011;
    
    
    koff = a.*exp(b.*f);

end

