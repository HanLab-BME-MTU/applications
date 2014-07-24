function w=FFTwaveletanalysis2D(im,FA1,FA2,J);

% FFTWAVELETANALYSIS2D FFT-based implementation of the wavelet transform.
% 	w=FFTwaveletanalysis2D(im,FFTanalysisfilters1,FFTanalysisfilters2,J) 
%	computes the wavelet transform of a volume using a Fourier method.  
%	It uses periodic boundary conditions.  The length of the signal must 
%	be a power of two and the frequency responses of the filters specified 
%	in FFTanalysisfilters.
% 	
% 	Input:
% 	im=input image, of size 2^N=length of FFTanalysisfilters
% 	J=depth of the decomposition, i.e., J wavelet bands (wav1 to wavJ) 
% 	+ 1 lowpass (lowJ) FFTanalysisfilters=[lowpassfilter;highpassfilter]
% 	
% 	See also FFTWAVELETSYNTHESIS, FFTFRACTSPLINEFILTERS, WEXTRACT.
% 	
% 	Authors: Thierry Blu, Dimitri Van De Ville
% 	Biomedical Imaging Group, EPFL, Lausanne, Switzerland.
% 	This software is downloadable at http://bigwww.epfl.ch/
% 	
% 	References:
% 	[1] M. Unser and T. Blu, "Fractional splines and wavelets," 
% 	SIAM Review, in press.
% 	[2] M. Unser and T. Blu, "Construction of fractional spline wavelet bases," 
% 	Proc. SPIE vol 3813, Wavelet Applications in Signal and Image 
% 	Processing VII, in press. 

FA.x=FA1; FA.y=FA2;

[M,N]=size(im);

if length(FA.x)~=M || length(FA.y)~=N
	disp(' ')
	disp('The size of the input signal and of the filters must match!')
	disp(' ')
	w=[];
	return
end

G1=FA.x(1,:); G2=FA.y(1,:); 
H1=FA.x(2,:); H2=FA.y(2,:); 

w=zeros(M,N);
for j=1:J,
	%
	% Computation of the outputs y and z
	%
	for i=1:M,
		X=fft(im(i,:),N);
		Y=G2.*X;
		Z=H2.*X;
		Y=1/2*(Y(1:N/2)+Y(N/2+(1:N/2)));
		Z=1/2*(Z(1:N/2)+Z(N/2+(1:N/2)));
		z=ifft(Z,N/2);
		im(i,:)=real([ifft(Y,N/2) z]);
		%tp1=real([ifft(Y,N/2) z]);
	end;
	for i=1:N,
		X=fft(reshape(im(:,i),1,M),M);
		Y=G1.*X;
		Z=H1.*X;
		Y=1/2*(Y(1:M/2)+Y(M/2+(1:M/2)));
		Z=1/2*(Z(1:M/2)+Z(M/2+(1:M/2)));
		z=ifft(Z,M/2);
		im(:,i)=real([ifft(Y,M/2) z]');
		%tp=real([ifft(Y,M/2) z]);
	end
	w(1:M,1:N)=im;
	M=M/2;N=N/2;
	im=w(1:M,1:N);
    
  
	G1=G1(1:2:length(G1)); G2=G2(1:2:length(G2)); 
	H1=H1(1:2:length(H1)); H2=H2(1:2:length(H2)); 
end;


