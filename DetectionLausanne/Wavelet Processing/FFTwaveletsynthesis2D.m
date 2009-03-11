function x=FFTwaveletsynthesis2D(w,FS1,FS2,J);

% FFTWAVELETSYNTHESIS2D
% 	
% 	Biomedical Imaging Group, EPFL, Lausanne, Switzerland.
% 	This software is downloadable at http://bigwww.epfl.ch/
% 	

[M,N]=size(w);

if length(FS1)~=M || length(FS2)~=N
	disp(' ')
	disp('The size of the input signal and of the filters must match!')
	disp(' ')
	w=[];
	return
end

%
% Reconstruction of the signal from its
% bandpass components
%

G1=conj(FS1(1,:)); G2=conj(FS2(1,:)); 
H1=conj(FS1(2,:)); H2=conj(FS2(2,:)); 

M=M/2^J;
N=N/2^J;

for j=J:-1:1,
	Hd1=H1(1:2^(j-1):length(H1));
	Hd2=H2(1:2^(j-1):length(H2));
	Gd1=G1(1:2^(j-1):length(G1));
	Gd2=G2(1:2^(j-1):length(G2));

	im=w(1:2*M,1:2*N);

  	for i=1:2*N,
	  Y=fft(reshape(im(1:M,i),1,M),M);
	  z=reshape(im((M+1):2*M,i),1,M);
	  Z=fft(z,M);
	  Y0=Gd1(1:M).*Y + Hd1(1:M).*Z;
	  Y1=Gd1(M+(1:M)).*Y + Hd1(M+(1:M)).*Z;
	  Y=[Y0 Y1];
	  im(:,i)=reshape(real(ifft(Y,2*M)),2*M,1);
	end	
	for i=1:2*M,
	  Y=fft(im(i,1:N),N);
	  z=im(i,(N+1):2*N);
	  Z=fft(z,N);
	  Y0=Gd2(1:N).*Y + Hd2(1:N).*Z;
	  Y1=Gd2(N+(1:N)).*Y + Hd2(N+(1:N)).*Z;
	  Y=[Y0 Y1];
	  im(i,:)=real(ifft(Y,2*N));    
	end

	M=2*M;
	N=2*N;

	w(1:M,1:N)=im;
  	x=w;
end

