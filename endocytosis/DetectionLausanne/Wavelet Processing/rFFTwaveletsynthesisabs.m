function x=rFFTwaveletsynthesisabs(w,FFTsynthesisfilters,J);

% FFTWAVELETANALYSIS FFT-based implementation of the inverse wavelet 
% 	transform.
% 	x=FFTwaveletsynthesis(w,FFTsynthesisfilters,J) computes the inverse 
% 	wavelet transform of w.  This function is the inverse of 
% 	FFTwaveletanalysis.  It uses periodic boundary conditions.  The 
% 	wavelet coefficient vector has the following fine-to-coarse 
% 	organization: w=[wav1 wav2 ...  wavJ lowJ]
% 
% 	Input:
% 	w=wavelet transform, of size 2^N=length of FFTsynthesisfilters
% 	FFTsynthesisfilters=[lowpassfilter;highpassfilter]
% 	J=depth of the decomposition, i.e., J wavelet bands + 1 lowpass
% 	
% 	Output:
% 	x=signal of size 2^N
% 	
% 	See also FFTWAVELETANALYSIS, FFTFRACTSPLINEFILTERS, WEXTRACT.
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

M0=length(w);
M=length(w)/(J+1);

if length(FFTsynthesisfilters)~=M
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

FFTsynthesisfilters=FFTsynthesisfilters/sqrt(2);

G=conj(FFTsynthesisfilters(1,:));
H=conj(FFTsynthesisfilters(2,:));

tmp=ones(1,M);
for i=J:-1:1,
    step=2^(i-1);
    tmp=tmp.*G(mod(0:step:M*step-1,M)+1);
end
tmp=convertabs(tmp);

LP=fft(w(M0-M+1:M0)).*tmp;
M0=M0-M;

for j=J:-1:1

	step=2^(j-1);

	tmp=H(mod(0:step:M*step-1,M)+1);

        for jtmp=j-1:-1:1,
	    step=2^(jtmp-1);
            tmp=tmp.*G(mod(0:step:M*step-1,M)+1);
            %since now we only need the lowpass filter  
        end
        tmp=convertabs(tmp);

	HP=fft(w(M0-M+1:M0)).*tmp;
	M0=M0-M;

	LP=LP+HP;

end

x=real(ifft(LP));

function v=convertabs(a)
v=fft(abs(ifft(a)));
