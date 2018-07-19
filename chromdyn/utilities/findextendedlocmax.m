function idx=findextendedlocmax(x)
%FINEXTENDEDLOCMAX returns index to true and quasi locmax
%
%

% c: 18/11/2002 dT

idx = [];
flag=1;
ti=1;
for i = 2:length(x)
   if(x(i)>x(i-1))
       flag=1;
       ti=i;
       %if last fill directly
       if i==length(x)
           idx=[idx,i];
       end;
   elseif flag & (x(i)<x(i-1) | i==length(x)) 
       flag=0;
       idx = [idx,ti];
       % fill both ends after a maximum horizontal line
       if (ti~=(i-1))
           idx=[idx,(i-1)];
       end;
   end;
end;