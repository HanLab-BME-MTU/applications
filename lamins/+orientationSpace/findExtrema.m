function [ out ] = findExtrema( orientationMatrix )
%orientationSpace.findExtrema finds the local maxima and minima by solving for the
%roots of the first derivative

    import orientationSpace.*;

    s = size(orientationMatrix);
    M = reshape(orientationMatrix,s(1)*s(2),s(3));
    % todo, deal with even case
    n = (s(3))/2;
%     n = n+0.5;

    
%     x = [0:n -n:-1];
    x = wraparoundN(0:s(3)-1,-n,n);
    % toeplitz?
    xx = wraparoundN(bsxfun(@minus,x,(0:(2*n-1))'),[-n n]);
    % xx = arrayfun(@(k) circshift(x,k,2),0:s(3),'UniformOutput',false);
    % xx = vertcat(xx{:});

    A = 2*exp(-xx.^2);
    [~,ind] = max(real(M),[],2);
%     xind = x(ind);

    M = real(M)/A*1e3;
%     out = zeros(s(1),s(2));
%     progressText(0,'Solving');
    nM = size(M,1);
    options = optimoptions('lsqnonlin','Jacobian','on','Display','off');
%     out = lsqnonlin(@orientationDerivWrapM,x(ind(:)),x(ind(:))-0.5,x(ind(:))+0.5,options);
%     out = reshape(out,s(1),s(2));
%     options = optimoptions('fsolve','Jacobian','on','Display','off');
    batchSize = s(1);
    nBatches = nM/batchSize;
%     Mcell = distributed.cell(1,nBatches);
%     xind = distributed.cell(1,nBatches);
% %     idx = distributed.zeros(nBatches,nBatches);
% %     out = cell(nM/batchSize,1);
%     spmd
%         for i=drange(1:nBatches)
%              idx = (1:batchSize)+(i-1)*batchSize;
%              Mcell(i) = {M(idx,:)};
%              xind(i) = {x(ind(idx))'};
%         end
%     end

    xind = x(ind)';
    


    Mcell = cell(1,nBatches);
    xindc = cell(1,nBatches);
    for i=1:nBatches
        idx = (1:batchSize)+(i-1)*batchSize;
        Mcell{i} = M(idx,:);
        xindc{i} = xind(idx);
    end
    
%     Mcell = distributed(Mcell);
%     xind = distributed(xind);

    f = getFunction(x,n,options);
    spmd
%         cod = codistributor1d(1);
%         M = codistributed(M,cod);
%         xind = codistributed(xind,cod);
% %         f = @(k) orientationDerivWrap(k,x,getLocalPart(M),n);
% %         offset = getLocalPart(xind);
% %         out = lsqnonlin(f,offset,offset-0.5,offset+0.5,options);
% %         size(getLocalPart(M))
% %         size(getLocalPart(xind))
%         Ml = getLocalPart(M);
%         xindl = getLocalPart(xind);
%         nParts = floor(size(Ml,1)/1024);
%         out = cell(nParts+1,1);
%         for i=1:nParts
%             idx = (1:1024)+(i-1)*1024;
%             out{i} = f(Ml(idx,:),xindl(idx));
%         end
%         out{end} = f(Ml((1+i*1024):end,:),xindl((1+i*1024):end,:));

        cod = codistributor1d(2);
        Mcell = codistributed(Mcell,cod);
        xindc = codistributed(xindc,cod);
        Ml = getLocalPart(Mcell);
        xl = getLocalPart(xindc);
        out = cell(size(Ml,2),1);
        for i=1:size(Ml,2)
            out{i} = f(Ml{i},xl{i});
        end
    end

%     parfor i=1:nM/batchSize
% %         disp(i);
% %         offset = xind(i);
%           offset = xindc{i};
% %         xxx = wraparoundN(x-offset,-n,n);
% %         f = @(o) -2*M(i,:)*exp(-(xxx-o).^2)';
% %         out(i) = fminbnd(f,-0.5,0.5) + x(ind(i));
% %         f = @(x) orientationDerivWrap(x,xxx,M(i,:));
%           f = @(k) orientationDerivWrap(k,x,Mcell{i},n);
% %         out(i) = fsolve(f,0,options)+offset;
% %         out(i) = fzero(f,0) + offset;
% %           out(i) = lsqnonlin(f,0,-0.5,0.5,options) + offset;
%           out{i} = lsqnonlin(f,offset,offset-0.5,offset+0.5,options);
% %         progressText(i/nM*1024);
%     end
%     doIt = getFunction(x,n,options);
%     out = cellfun(doIt,Mcell,xind,'UniformOutput',false);
%     out = gather(out);
%     out = horzcat(out{:});

%     function [F,J] = orientationDerivWrapM(k)
%         [F,J] = orientationDerivWrap(k,x,M,n);
%     end
%     function out = doIt(M,offset)
%         out = lsqnonlin(@doIt2,offset,offset-0.5,offset+0.5,options);
%         function [F,J] = doIt2(k)
%             [F,J] = orientationDerivWrap(k,x,M,n);
%         end
%     end

end
function f = getFunction(x,n,options)
    f = @doIt;
    function out = doIt(M,offset)
        out = lsqnonlin(@doIt2,offset,offset-0.5,offset+0.5,options);
        function [F,J] = doIt2(k)
            [F,J] = orientationDerivWrap(k,x,M,n);
        end
    end
end
function [F,J] = orientationDerivWrap(k,x,M,n)
    t = wraparoundN(bsxfun(@minus,k,x),-n,n);
    et = exp(-t.^2);
    F = sum(M.*(t.*et),2);
    if(nargout > 1)
        J = sum(M.*(et -2*t.^2.*et),2);
        J = spdiags(J(:),0,length(F),length(k));
    end
end
function [F,J] = orientationDeriv(x,n,a)
    t = x - n;
    et = exp(-t.^2);
    F = sum(a.*t.*et);
    if(nargout > 1)
        J = sum(a.*(et -2*t.^2.*et));
    end
end
function W = jmfunc(Jinfo,Y,flag)
    if(~flag)
        Jinfo = Jinfo.^2;
    end
    W = bsxfun(@times,Jinfo(:),Y);
end
