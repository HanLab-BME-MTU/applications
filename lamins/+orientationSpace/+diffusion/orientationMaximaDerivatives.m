function [ dnm_dKn ] = orientationMaximaDerivatives( rho, K, derivOrder, lm )
%ORIENTATIONMAXIMADERIVATIVES Find the K derivatives of local maxima
%
% INPUT
% rho - regularly spaced samples of orientation response at K,
%       nSamples x numel(K)
% lm  - orientation local maxima of rho, via interpft_extrema
%       maxNLocalMaxima x numel(K)
% K - angular order, vector of K values correpsond to rho and lm
% derivOrder - highest derivative desired
%
% OUTPUT
% derivatives

    D = 2*pi^2;

%% Calculate derivative in terms of t
%t = 1./(2*K+1).^2;

% n = derivOrder
% m = theta_{m} (local maxima location)

if(nargin < 4)
    lm = interpft_extrema(rho);
    lm = orientationSpace.diffusion.alignExtrema(lm);
end

rho_derivs = interpft1_derivatives(rho,lm,2:derivOrder*2+1);

[dnm_dtn,maximaDerivatives] = dqm_dtq(derivOrder,D,rho_derivs);
    
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

for d = 1:derivOrder
    dnm_dKn(:,:,d) = translate_from_t_to_K(d,cat(3,maximaDerivatives{:}),K);
end


end

function deriv = total_dq_dtq_partial_dnrho_dmn(q,n,D,rho_derivs,maximaDerivatives)
    assert(q >= 0);
    assert(n > 0);
    if(n == 1)
        % Total derivative of 1st partial derivative with respect to
        % orientation. Always zero by definition of orientation local
        % maxima.
        deriv = 0;
%         fprintf('GET   total q=%d, n=%d\n',q,n);
%         fprintf('END total q=%d, n=%d\n',q,n);
        return;
    end
    if(q == 0)
        % No total derivative, answer is just the nth partial derivative with
        % respect to orientation
        deriv = rho_derivs(:,:,n-1);
%         fprintf('END total q=%d, n=%d\n',q,n);
%         fprintf('GET   total q=%d, n=%d\n',q,n);
        return;
    end
%     fprintf('START total q=%d, n=%d\n',q,n);

    deriv = 0;
    for l = 1:q
        binom = nchoosek(q-1,l-1);
%         deriv = deriv + binom.*total_dq_dtq_partial_dnrho_dmn(q-l,n+1,D,rho_derivs).*dqm_dtq(l,D,rho_derivs);
        deriv = deriv + binom.*total_dq_dtq_partial_dnrho_dmn(q-l,n+1,D,rho_derivs,maximaDerivatives).*maximaDerivatives{l};
    end
    deriv = deriv +  D * total_dq_dtq_partial_dnrho_dmn(q-1,n+2,D,rho_derivs,maximaDerivatives);
%     fprintf('END   total q=%d, n=%d\n',q,n);
end

function [deriv,maximaDerivatives] = dqm_dtq(q,D,rho_derivs,maximaDerivatives)
    % order of the partial derivative is 1
    % need rho_derivs up to 1+q*2;
%     fprintf('START dqm_dtq q=%d\n',q);
    assert(q > 0);
    n = 1;
    deriv = 0;
    if(nargin < 4)
        maximaDerivatives = cell(1,q);
        for l=1:q-1
            maximaDerivatives{l} = dqm_dtq(l,D,rho_derivs,maximaDerivatives);
        end
    end
    for l=1:q-1
        binom = nchoosek(q-1,l-1);
        deriv = deriv + binom.*total_dq_dtq_partial_dnrho_dmn(q-l,n+1,D,rho_derivs,maximaDerivatives).*maximaDerivatives{l};
    end
    deriv = deriv + D * total_dq_dtq_partial_dnrho_dmn(q-1,n+2,D,rho_derivs,maximaDerivatives);
    % Divide by second partial derivative with respect to orientation
    deriv = -deriv./rho_derivs(:,:,1);
    maximaDerivatives{q} = deriv;
%     fprintf('END dqm_dtq q=%d\n',q);
end

function tpm = total_partial_matrix()
end

function dqm_dKq = translate_from_t_to_K(q,dqm_dtq_v,K)
    %% Calculate Faa di Bruno coefficients
    part = partitions(q);
    faa_di_bruno = factorial(q);
    faa_di_bruno = faa_di_bruno./prod(bsxfun(@power,factorial(1:q),part).');
    faa_di_bruno = faa_di_bruno./prod(factorial(part).');
    %% Calculate order of the local maxima derivative
    derivOrder = sum(part.');
    %% Calculate derivatives of t with respect to K
    n = shiftdim(1:q,-1);
    dnt_dKn = factorial(n+1) .* (-2).^(n );
    dnt_dKn = bsxfun(@times,dnt_dKn, ...
        bsxfun(@power,1./(2*K+1),n+2));
    dnt_dKn_pow = bsxfun(@power,dnt_dKn,shiftdim(part.',-2));
    dnt_dKn_pow = prod(dnt_dKn_pow,3);
    dqm_dKq = bsxfun(@times,dnt_dKn_pow,shiftdim(faa_di_bruno,-2));
    dqm_dKq = reshape(dqm_dKq,[1 size(dqm_dKq,2) size(dqm_dKq,4)]);
    dqm_dKq = bsxfun(@times,dqm_dKq,dqm_dtq_v(:,:,derivOrder));
    dqm_dKq = sum(dqm_dKq,3);
end