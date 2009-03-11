function x=FFTwaveletsynthesis2Dabs(w,FS1,FS2,J,PA);

% FFTWAVELETSYNTHESIS2DABS
%       FFT-based implementation of the inverse 2D wavelet transform with
%       absolute value wavelet
%

[Mx,My]=size(w);

if length(FS1)~=Mx || length(FS2)~=My
	disp(' ')
	disp('The size of the input signal and of the filters must match!')
	disp(' ')
	w=[];
	return
end

PS=1./PA;

Hx=conj(FS1(1,:)); Hy=conj(FS2(1,:));
Gx=conj(FS1(2,:)); Gy=conj(FS2(2,:));

sz1=Mx;
sz2=My;

Mx=Mx/2;
My=My/2;

mx=sz1;
my=sz2;

x=zeros(size(w));

Hcx=ones(size(Hx));
Hcy=ones(size(Hy));

for j=1:J,

	if j==1,
	  G1x=Gx;
	  G1y=Gy;
	else
	  G1x(1:sz1/2)=G1x(1:2:sz1);
	  G1x(sz1/2+1:sz1)=G1x(1:sz1/2);

	  G1y(1:sz2/2)=G1y(1:2:sz2);
	  G1y(sz2/2+1:sz2)=G1y(1:sz2/2);
	end;
	
	if j==1,
	  H1x=ones(size(Hx));
	  H1y=ones(size(Hy));

	  H2x=Hx;
	  H2y=Hy;
	else
	  H1x=H2x;
	  H1y=H2y;

	  H2x(1:sz1/2)=H2x(1:2:sz1);
	  H2x(sz1/2+1:sz1)=H2x(1:sz1/2);

	  H2y(1:sz2/2)=H2y(1:2:sz2);
	  H2y(sz2/2+1:sz2)=H2y(1:sz2/2);
	end;

	Hcx=Hcx.*H1x;
	Hcy=Hcy.*H1y;

	% 1 HL
	z=w(mx-Mx+1:mx,1:My);
	Z=fftn(z); Z=duplicate3D(Z,j);
	F=construct3D(G1x.*Hcx,H2y.*Hcy,PS);
	RZ=Z.*F;
	x=x+real(ifftn(RZ));

	% 2 HH
	z=w(mx-Mx+1:mx,my-My+1:my);
	Z=fftn(z); Z=duplicate3D(Z,j);
	F=construct3D(G1x.*Hcx,G1y.*Hcy,PS);
	RZ=Z.*F;
	x=x+real(ifftn(RZ));

	% 4 LH
	z=w(1:Mx,my-My+1:my);
	Z=fftn(z); Z=duplicate3D(Z,j);
	F=construct3D(H2x.*Hcx,G1y.*Hcy,PS);
	RZ=Z.*F;
	x=x+real(ifftn(RZ));

	mx=mx-Mx; my=my-Mx; 
	Mx=Mx/2;  My=My/2;  
end

Mx=2*Mx; My=2*My;

z=w(1:Mx,1:My);
Z=fftn(z);

Z=duplicate3D(Z,j);

if j==1,
  H1x=Hx;
  H1y=Hy;
else
  H1x(1:sz1/2)=H1x(1:2:sz1);
  H1x(sz1/2+1:sz1)=H1x(1:sz1/2);

  H1y(1:sz2/2)=H1y(1:2:sz2);
  H1y(sz2/2+1:sz2)=H1y(1:sz2/2);
end;

Fx=H1x.*Hcx;
Fy=H1y.*Hcy;

F=construct3D(Fx,Fy,PS);
RZ=Z.*F;
x=x+real(ifftn(RZ));
return;

function Z=duplicate3D(Z,j);
for iter=1:j, 
  z1=size(Z,1); z2=size(Z,2);
  Z(z1+1:2*z1,:)=Z(1:z1,:);
  Z(:,z2+1:2*z2)=Z(:,1:z2);
end;
return

function F=construct3D(Fx,Fy,PS)
Fx=ifft(Fx); Fy=ifft(Fy); 
sz1=length(Fx); sz2=length(Fy);
F=Fx'*Fy;
F=reshape(F,sz1,sz2);
if size(PS,1)>0,      % postfilter? (compensation for prefilter)
  F=fftn(F);
  F=F.*PS;
  F=ifftn(F);
end;
F=fftn(abs(F)); % put absolute values of filter coefficients!
return
