function koff = koffcb(f,ion)
% Calculates koff as a function of force for a5b1 integrins (kong et al.
% JCB). Force in N, koff in s-1.

f = abs(f);

a = 0.*f;
b = 0.*f;

switch ion
    
    case 'mg'
    
    a(f <= 13.2e-12) = 0.1905;
    b(f <= 13.2e-12) = 2.333e+011;
    
    a(f > 30e-12) = 0.04527;
    b(f > 30e-12) = 7.251e+010;
    
    a(a ==0) = 30.94;
    b(b ==0) = -1.489e+011;
    
    koff = a.*exp(b.*f);

    case 'mg_earlyslipmore' % for arp2/3 inhibition
    
    a(f <= 0.5e-13) = 0.1905;
    b(f <= 0.5e-13) = 2.333e+010;
    
    a(f > 0.5e-13) = 0.1905;
    b(f > 0.5e-13) = 7e+011;
    
    a(a ==0) = 30.94;
    b(b ==0) = -1.489e+011;
    
    koff = a.*exp(b.*f);

    case 'mg_earlyslip' %for formin-inhibition
    
    a(f <= 0.3e-13) = 0.1905;
    b(f <= 0.3e-13) = 2.333e+011;
    
    a(f > 0.3e-13) = 0.1905;
    b(f > 0.3e-13) = 5.95e+011;
    
    a(a ==0) = 30.94;
    b(b ==0) = -1.489e+011;
    
    koff = a.*exp(b.*f);

    case 'mn'
    
    %     a(f <= 32.8e-12) = 6.237;
    %     b(f <= 32.8e-12) = -0.1466e12;
    %
    %     a(f > 32.8e-12) = 0.04137;
    %     b(f > 32.8e-12) = 5.262e+010;
    
    a = 11.27;       % 16.67; %
    b = -1.601e+011; %-1.945e11;
    c = 0.0008882;   % 0.001442;
    d = 1.225e+011;  % 1.356e11;
    koff = a.*exp(b.*f) + c.*exp(d.*f);
    
    case 'cm' %Custom curve to fit data
    
    c = 0.0008882;   % 0.001442;
    d = 1.225e+011;  % 1.356e11;
    a = 11.27;       % 16.67; %
    b = -1.601e+011; %-1.945e11;
    g = 1000; % Increase in koff at zero force wrt mn ion.
    fend = 0.15e-12; % Force at which the added curve vanishes.
    h = -log(0.1./g)./fend;
    
    % b = -1.75e11;
    % a = -(c.*d./b)*exp((d-b).*34e-12);
    koff = a.*exp(b.*f) + c.*exp(d.*f) + g.*exp(-h.*f);
    
    case 'b3' % Data for b3 integrin
        
    a =      0.4012; 
    b = -3.488e+010; 
    c =   0.0004173;  
    d =  1.824e+011;
    
    koff = a.*exp(b.*f) + c.*exp(d.*f);
    
    case 'b3_lateslip' % Data for b3 integrin
        
    a =      0.8; %0.4012; 
    b = -3.488e+010; 
    c =   0.0004173;  
    d =  1.824e+10;
    
    koff = a.*exp(b.*f) + c.*exp(d.*f);

    case 'b3cm' % Data for b3 integrin, custom
        
    a =      0.4012; 
    b = -3.488e+010; 
    c =   0.0004173;  
    d =  1.824e+011;
        g = 300; % Increase in koff at zero force wrt mn ion.
    fend = 0.15e-12; % Force at which the added curve vanishes.
    h = -log(0.1./g)./fend;
    
    koff = a.*exp(b.*f) + c.*exp(d.*f) + g.*exp(-h.*f);
end

end

