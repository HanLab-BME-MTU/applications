function A2=showdwt2d(A,J,varargin);

m1=size(A,1);
n1=m1/2+1;
m2=size(A,2);
n2=m2/2+1;

A2=A;

for iter=1:J,
    if length(varargin)==0,
    A2(n1:m1,n2:m2)=rescale(A(n1:m1,n2:m2));
    A2(1:n1-1,n2:m2)=rescale(A(1:n1-1,n2:m2));
    A2(n1:m1,1:n2-1)=rescale(A(n1:m1,1:n2-1));
    end;
    A2(n1,1:m2)=1.5;
    A2(1:m1,n2)=1.5;
    m1=n1-1;
    n1=(n1-1)/2+1;
    m2=n2-1;
    n2=(n2-1)/2+1;
end;

A2(1:m1,1:m2)=rescale(A(1:m1,1:m2));

function R=rescale(sub);
%R=(sub-min(sub(:)))/abs(max(sub(:))-min(sub(:)));
R=sub/max(abs(sub(:)));
