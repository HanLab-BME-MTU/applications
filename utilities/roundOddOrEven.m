function out=roundOddOrEven(x,oddOrEven,infOrZero)
%rounds to the next even or odd number
%
%SYNOPSIS out=roundOddOrEven(x,oddOrEven,infOrZero)
%
%INPUT    x: vector or matrix of values
%         optional:
%         oddOrEven: ['odd'],'even': if round to odd or even number [default]
%         infOrZero: 'inf',['close'],'zero': if round towards infinity, closest
%         number or zero [default]
%
%OUTPUT   out: rounded values
%
%created: 2/03, Jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%test input
switch nargin
    case 1
        oddOrEven='odd';
        infOrZero='';
    case 2
        infOrZero='';
        if ~strcmp(oddOrEven,'odd')&~strcmp(oddOrEven,'even')
            error('wrong input for oddOrEven');
            return
        end
    case 3
        if ~strcmp(oddOrEven,'odd')&~strcmp(oddOrEven,'even')
            error('wrong input for oddOrEven');
            return
        end
        if ~strcmp(infOrZero,'inf')&~strcmp(infOrZero,'zero')&~strcmp(infOrZero,'close')
            error('wrong input for infOrZero');
            return
        end
end

%evaluation switch
switch infOrZero
    case 'inf'
        switch oddOrEven
            case 'odd'
                %remember sign
                sig=sign(x);
                %round
                x=ceil(sig.*x).*sig;
                %make odd
                out=x+(1-mod(x,2)).*sig;
            case 'even'
                %remember sign
                sig=sign(x);
                %round
                x=ceil(sig.*x).*sig;
                %make even
                out=x+(mod(x,2)).*sig;
        end
    case 'zero'
        switch oddOrEven
            case 'odd'
                %remember sign
                sig=sign(x);
                %round
                x=fix(x);
                %make odd
                %first: check for zeros
                zeroIdx=find(x==0);
                if ~isempty(zeroIdx)
                    x(zeroIdx)=sig(zeroIdx);
                end
                out=x-(1-mod(x,2)).*sig;
            case 'even'
                %remember sign
                sig=sign(x);
                %round
                x=fix(x);
                %make even
                out=x-(mod(x,2)).*sig;
        end
    case 'close'
        switch oddOrEven
            case 'odd'
                %remember sign
                sig=sign(x);
                %remember modulo
                modulo=floor(mod(x.*sig,2));
                %round
                x=fix(x);
                %make odd
                out=x+(1-modulo).*sig;
            case 'even'
                %remember sign
                sig=sign(x);
                %remember modulo
                modulo=floor(mod(x.*sig,2));
                %round
                x=fix(x);
                %make even
                out=x+(modulo).*sig;
        end
end