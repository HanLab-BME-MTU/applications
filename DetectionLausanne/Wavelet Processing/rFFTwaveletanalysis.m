function w=rFFTwaveletanalysis(x,FFTanalysisfilters,J);

% FFTWAVELETANALYSIS FFT-based implementation of the wavelet transform.
% 	w=FFTwaveletanalysis(x,FFTanalysisfilters,J) computes the wavelet 
% 	transform of a signal x using a Fourier method.  It uses periodic 
% 	boundary conditions.  The length of the signal must be a power of 
% 	two and the frequency responses of the filters specified in FFTanalysisfilters.
% 	
% 	Input:
% 	x=input signal, of size 2^N=length of FFTanalysisfilters
% 	J=depth of the decomposition, i.e., J wavelet bands (wav1 to wavJ) 
% 	+ 1 lowpass (lowJ) FFTanalysisfilters=[lowpassfilter;highpassfilter]
% 	
% 	Output:
% 	w=[wav1 wav2 ...  wavJ lowJ]: vector of pooled wavelet 
% 	coefficients, size 2^N
% 	
% 	See also FFTWAVELETSYNTHESIS, FFTFRACTSPLINEFILTERS, WEXTRACT.
% 	
% 	Author: Thierry Blu, October 1999
% 	Biomedical Imaging Group, EPFL, Lausanne, Switzerland.
% 	This software is downloadable at http://bigwww.epfl.ch/
% 	
% 	References:
% 	[1] M. Unser and T. Blu, "Fractional splines and wavelets," 
% 	SIAM Review, Vol. 42, No. 1, pp. 43--67, January 2000.
% 	[2] M. Unser and T. Blu, "Construction of fractional spline wavelet bases," 
% 	Proc. SPIE, Wavelet Applications in Signal and Image Processing VII,
%     Denver, CO, USA, 19-23 July, 1999, vol. 3813, pp. 422-431. 
% 	[3] T. Blu and M. Unser, "The fractional spline wavelet transform: definition and 
%	implementation," Proc. IEEE International Conference on Acoustics, Speech, and 
%	Signal Processing (ICASSP'2000), Istanbul, Turkey, 5-9 June 2000, vol. I, pp. 512-515 .

M=length(x);

if length(FFTanalysisfilters)~=M
	disp(' ')
	disp('The size of the input signal and of the filters must match!')
	disp(' ')
	w=[];
	return
end

FFTanalysisfilters=FFTanalysisfilters/sqrt(2);

% Fourier transform of the signal
LP=fft(x);

G=FFTanalysisfilters(1,:);
H=FFTanalysisfilters(2,:);

w=[];

for j=1:J

	step=2^(j-1);

	X=LP.*G(mod(0:step:M*step-1,M)+1);
	HP=LP.*H(mod(0:step:M*step-1,M)+1);

	z=real(ifft(HP));
	w=[w z];
	
	LP=X;
end

w=real([w ifft(LP)]);

