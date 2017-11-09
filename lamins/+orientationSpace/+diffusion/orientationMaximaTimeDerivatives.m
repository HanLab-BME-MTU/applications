function [ dtn_dnm, maximaDerivatives, dnK_dmn ] = orientationMaximaTimeDerivatives( rho, K, derivOrder, lm , period, freq)
%orientationMaximaTimeDerivatives Find the derivative of time/K with
%respect to orientation of local orientaiton maxima
%
% INPUT
% rho - regularly spaced samples of orientation response at K,
%       nSamples x numel(K)
% K - angular order, vector of K values correpsond to rho and lm
%     If scalar, then K is expanded to match size(rho,2)
% derivOrder - highest derivative desired
% lm  - orientation local maxima of rho, via interpft_extrema
%       maxNLocalMaxima x numel(K)
%
% OUTPUT
% full derivatives of the orientation local maxima with respect to K

% Mark Kittisopikul, August 2017

    if(nargin < 5 || isempty(period))
        period = 2*pi;
    end
    if(nargin < 6 || isempty(freq))
        freq = false;
    end
    
    D = period.^2/2;
    
    % For period of 2*pi
%     D = 2*pi^2;
    % For period of pi
    % D = pi^2/2;

%% Calculate derivative in terms of t
%t = 1./(2*K+1).^2;

% n = derivOrder
% m = theta_{m} (local maxima location)

if(nargin < 4 || isempty(lm))
    lm = interpft_extrema(rho);
    lm = orientationSpace.diffusion.alignExtrema(lm,period);
end

derivOrders = 2:derivOrder*2+1;

if(size(rho,3) == length(derivOrders))
    % rho includes derivative information
    method = 'horner';
    if(freq)
        method = 'horner_freq';
    end
    rho_derivs = interpft1([0 period],rho,lm,method,10,false);
else
    if(~ismatrix(rho))
    %     rho_sz = size(rho);
        rho = rho(:,:);
    end
    if(~ismatrix(lm))
        lm_sz = size(lm);
        lm = lm(:,:);
    end

    rho_derivs = interpft1_derivatives(rho,lm,2:derivOrder*2+1,period,freq);
    
end

[maximaDerivatives] = dqt_dtm(derivOrder,D,rho_derivs);

dtn_dnm = maximaDerivatives{end};
    
% if(derivOrder == 1)
% %     rho_derivs = interpft1_derivatives(rho,lm,[2 3]);
%     dnm_dtn = -D.*rho_derivs(:,:,2)./rho_derivs(:,:,1);
% else
% end

%% Calculate derivatives of t by K

% n = shiftdim(1:derivOrder,-1);
% dnt_dKn = factorial(n+1) .* (-2).^(n );
% dnt_dKn = bsxfun(@times,dnt_dKn, ...
%     bsxfun(@power,1./(2*K+1),n+2));
% keyboard

%% Translate derivative with respect to t to with respect to K

% if(derivOrder == 1)
% %     keyboard;
%     dnm_dKn = bsxfun(@times,dnm_dtn,dnt_dKn(:,:,1));
% else
%     dnm_dKn = dnm_dtn;
% end

if(nargout > 2)

    dnK_dtn = get_K_derivatives_with_respect_to_t(derivOrder,K);
    K_expansion_factor = 1;
    if(isscalar(K))
        K_expansion_factor = size(rho,2);
    end
    dnK_dtn = repmat(dnK_dtn,[size(lm,1) K_expansion_factor 1]);

    dnK_dmn = zeros([size(lm),derivOrder]);
    for d = 1:derivOrder
        dnK_dmn(:,:,d) = translate_from_t_to_K_hard(d,cat(3,maximaDerivatives{:}),K,dnK_dtn);
    %     dnm_dKn(:,:,d) = translate_from_t_to_K(d,cat(3,maximaDerivatives{:}),K);
    end

    if(exist('lm_sz','var'))
        dnK_dmn = reshape(dnK_dmn,[lm_sz derivOrder]);
    end

end



end

function [maximaDerivatives] = dqt_dtm(q,D,rho_derivs,maximaDerivatives)

if(nargin < 4)
    maximaDerivatives = cell(1,q);
    if(q > 1)
        maximaDerivatives(1:q-1) = dqt_dtm(q-1,D,rho_derivs);
    end
end

switch(q)
    case 1
        maximaDerivatives{q} =       rho_derivs(:,:,2-1);
    case 2
        maximaDerivatives{q} =       rho_derivs(:,:,3-1) ...
                             + 2*D*  rho_derivs(:,:,4-1).*maximaDerivatives{1} ...
                             + D^2*  rho_derivs(:,:,5-1).*(maximaDerivatives{1}).^2;
    case 3
        maximaDerivatives{q} =       rho_derivs(:,:,4-1).*(1+3*D*maximaDerivatives{2}) ...
                             + 3*D*  rho_derivs(:,:,5-1).*maximaDerivatives{1} ...
                             + 3*D^2*rho_derivs(:,:,5-1).*maximaDerivatives{1}.*maximaDerivatives{2} ...
                             + 3*D^2*rho_derivs(:,:,6-1).*(maximaDerivatives{1}).^2 ...
                             + D^3*  rho_derivs(:,:,7-1).*(maximaDerivatives{1}).^3;
    case 4
        maximaDerivatives{q} = 4*D*  rho_derivs(:,:,4-1).*maximaDerivatives{3} ...
                             + 4*D^2*rho_derivs(:,:,5-1).*maximaDerivatives{1}.*maximaDerivatives{3} ...
                             +       rho_derivs(:,:,5-1) ...
                             + 6*D*  rho_derivs(:,:,5-1).*maximaDerivatives{2} ...
                             + 3*D^2*rho_derivs(:,:,5-1).*(maximaDerivatives{2}).^2 ...
                             + 4*D*  rho_derivs(:,:,6-1).*maximaDerivatives{1} ...
                             +12*D^2*rho_derivs(:,:,6-1).*maximaDerivatives{1}.*maximaDerivatives{2} ...
                             + 6*D^2*rho_derivs(:,:,7-1).*(maximaDerivatives{1}).^2 ...
                             + 6*D^3*rho_derivs(:,:,7-1).*(maximaDerivatives{1}).^2.*maximaDerivatives{2} ...
                             + 4*D^3*rho_derivs(:,:,8-1).*(maximaDerivatives{1}).^3 ...
                             + D^4*  rho_derivs(:,:,9-1).*(maximaDerivatives{1}).^4;
    otherwise
        import orientationSpace.diffusion.*;
        s = calcThetaDeriv(q,1);
        s = s(arrayfun(@(s) s.timed(1,q) == 0,s));
        ss = combineTerms(s);
        ss.timed = ss.timed(:,any(ss.timed));
        % Find partial derivative number
        [~,d] = find(ss.rhod);
        d = d - 1;
        front = ss.coeff.*D.^(ss.D);
        maximaDerivatives{q} = bsxfun(@times,shiftdim(front,-2),rho_derivs(:,:,d));
        for qq = 1:size(ss.timed,2)
            maximaDerivatives{q} = maximaDerivatives{q}.*bsxfun(@power,maximaDerivatives{qq},shiftdim(ss.timed(:,qq),-2));
        end
        maximaDerivatives{q} = sum(maximaDerivatives{q},3);
%         error('Not implemented');
        % See calcThetaDeriv
end
maximaDerivatives{q} = -maximaDerivatives{q}./rho_derivs(:,:,3-1)./D;


end

function dqK_dmq = translate_from_t_to_K_hard(q,dqt_dmq_v,K,dnK_dtn)
% Translate derivatives with respect to t to derivatives with respect to K
% Hardcoded version for efficiency
% q - scalar, order of the derivative with respect to K
% dqK_dmq - derivatives of K with respect to m
% K - 
    if(q > 6)
        dqK_dmq = translate_from_t_to_K(q,dnK_dtn,K);
        return;
    end
    
    dqt_dmq_v = repmat(dqt_dmq_v,[1 size(dqt_dmq_v,2)./size(dnK_dtn,2) 1]);
    
    switch(q)
        case 1
            % partitions(1)
            % 
            % ans =
            % 
            %      1
            dqK_dmq = bsxfun(@times,dnK_dtn(:,:,1),dqt_dmq_v(:,:,1));
        case 2
            % partitions(2)
            % 
            % ans =
            % 
            %      2     0
            %      0     1
            dqK_dmq = dnK_dtn(:,:,2).* dqt_dmq_v(:,:,1).^2 ... 
                    + dnK_dtn(:,:,1).* dqt_dmq_v(:,:,2);
        case 3
            % partitions(3)
            % 
            % ans =
            % 
            %      3     0     0
            %      1     1     0
            %      0     0     1
            dqK_dmq =   dnK_dtn(:,:,3) .* dqt_dmq_v(:,:,1).^3 ...
                    + 3*dnK_dtn(:,:,2) .* dqt_dmq_v(:,:,1)     .* dqt_dmq_v(:,:,2) ...
                    +   dnK_dtn(:,:,1) .* dqt_dmq_v(:,:,3);
        case 4
            % partitions(4)
            % 
            % ans =
            % 
            %      4     0     0     0
            %      2     1     0     0
            %      0     2     0     0
            %      1     0     1     0
            %      0     0     0     1
            dqK_dmq =   dnK_dtn(:,:,4) .* dqt_dmq_v(:,:,1).^4 ...
                    + 6*dnK_dtn(:,:,3) .* dqt_dmq_v(:,:,1).^2 .* dqt_dmq_v(:,:,2) ...
                    + 3*dnK_dtn(:,:,2) .* dqt_dmq_v(:,:,2).^2 ...
                    + 4*dnK_dtn(:,:,2) .* dqt_dmq_v(:,:,1)    .* dqt_dmq_v(:,:,3) ...
                    +   dnK_dtn(:,:,1) .* dqt_dmq_v(:,:,4);
        case 5
            % partitions(5)
            % 
            % ans =
            % 
            %      5     0     0     0     0
            %      3     1     0     0     0
            %      1     2     0     0     0
            %      2     0     1     0     0
            %      0     1     1     0     0
            %      1     0     0     1     0
            %      0     0     0     0     1
            dqK_dmq =   dnK_dtn(:,:,5) .* dqt_dmq_v(:,:,1).^5 ...
                    +10*dnK_dtn(:,:,4) .* dqt_dmq_v(:,:,1).^3 .* dqt_dmq_v(:,:,2) ...
                    +15*dnK_dtn(:,:,3) .* dqt_dmq_v(:,:,1)    .* dqt_dmq_v(:,:,2).^2 ...
                    +10*dnK_dtn(:,:,3) .* dqt_dmq_v(:,:,1).^2 .* dqt_dmq_v(:,:,3) ...
                    +10*dnK_dtn(:,:,2) .* dqt_dmq_v(:,:,2)    .* dqt_dmq_v(:,:,3) ...
                    + 5*dnK_dtn(:,:,2) .* dqt_dmq_v(:,:,1)    .* dqt_dmq_v(:,:,4) ...
                    +   dnK_dtn(:,:,1) .* dqt_dmq_v(:,:,5);
        case 6
            % partitions(6)
            % 
            % ans =
            % 
            %      6     0     0     0     0     0
            %      4     1     0     0     0     0
            %      2     2     0     0     0     0
            %      0     3     0     0     0     0
            %      3     0     1     0     0     0
            %      1     1     1     0     0     0
            %      0     0     2     0     0     0
            %      2     0     0     1     0     0
            %      0     1     0     1     0     0
            %      1     0     0     0     1     0
            %      0     0     0     0     0     1
            dqK_dmq =   dnK_dtn(:,:,6) .* dqt_dmq_v(:,:,1).^6 ...
                    +15*dnK_dtn(:,:,5) .* dqt_dmq_v(:,:,1).^4 .* dqt_dmq_v(:,:,2) ...
                    +45*dnK_dtn(:,:,4) .* dqt_dmq_v(:,:,1).^2 .* dqt_dmq_v(:,:,2).^2 ...
                    +15*dnK_dtn(:,:,3) .* dqt_dmq_v(:,:,2).^3  ...
                    +20*dnK_dtn(:,:,4) .* dqt_dmq_v(:,:,1).^3 .* dqt_dmq_v(:,:,3) ...
                    +60*dnK_dtn(:,:,3) .* prod(dqt_dmq_v(:,:,1:3),3) ...
                    +10*dnK_dtn(:,:,2) .* dqt_dmq_v(:,:,3).^2 ...
                    +15*dnK_dtn(:,:,3) .* dqt_dmq_v(:,:,1).^2 .* dqt_dmq_v(:,:,4) ...
                    +15*dnK_dtn(:,:,2) .* dqt_dmq_v(:,:,2) .* dqt_dmq_v(:,:,4) ...
                    + 6*dnK_dtn(:,:,2) .* dqt_dmq_v(:,:,1) .* dqt_dmq_v(:,:,5) ...
                    +   dnK_dtn(:,:,1) .* dqt_dmq_v(:,:,6);
        otherwise
            error('Should not have gotten here');
    end
end

function dqK_dmq = translate_from_t_to_K(q,dqt_dmq_v,K)
    %% Calculate Faa di Bruno coefficients
    part = partitions(q);
    faa_di_bruno = factorial(q);
    faa_di_bruno = faa_di_bruno./prod(bsxfun(@power,factorial(1:q),part).');
    faa_di_bruno = faa_di_bruno./prod(factorial(part).');
    %% Calculate order of the local maxima derivative
    derivOrder = sum(part.');
    %% Calculate derivatives of t with respect to K
    n = shiftdim(1:q,-1);
%     dnt_dKn = factorial(n+1) .* (-2).^(n );
%     dnt_dKn = bsxfun(@times,dnt_dKn, ...
%         bsxfun(@power,1./(2*K+1),n+2));
    dnK_dtn = get_K_derivatives_with_respect_to_t(n,K);
    % TODO, fix this, shaping 2017/08/03
    dqt_dmq_v_pow = bsxfun(@power,dqt_dmq_v,shiftdim(part.',-2));
    dqt_dmq_v_pow = prod(dqt_dmq_v_pow,3);
    dqK_dmq = bsxfun(@times,dqt_dmq_v_pow,shiftdim(faa_di_bruno,-2));
    
    dqK_dmq = reshape(dqK_dmq,[1 size(dqK_dmq,2) size(dqK_dmq,4)]);
    dqK_dmq = bsxfun(@times,dqK_dmq,dnK_dtn(:,:,derivOrder));
    dqK_dmq = sum(dqK_dmq,3);
end

function dnK_dtn = get_K_derivatives_with_respect_to_t(q,K)
        dnK_dtn = zeros(1,length(K),q);
        dnK_dtn(:,:,1) = (2*K+1).^3./-4;
        for i=2:q
            dnK_dtn(:,:,i) = dnK_dtn(:,:,i-1).*(2*K+1).^2/-2*(i*2-1);
        end
end

