function [frequencyAxis,singleSidedSpectrum] = accelerometerMeasure(time, rate);

% A program used to measure vibrations with an accelerometer.  The original
% .m file was lost during the reconfiguration of a computer.  
% We use the NI cDAQ-9171 to communicate with the NI 9234, cDAQ1Mod1.
% Also, you must run the accelerometer a few times to get a reasonable
% answer.  Presumably due to built up capacitance in the piezoelectric
% crystal.
% Mainly based off of FFT example from Matlab: http://www.mathworks.com/help/matlab/ref/fft.html
% Kevin Dean, 2015-12-02.
% 

daq.getDevices;

s = daq.createSession('ni');

addAnalogInputChannel(s,'cDAQ1Mod1',0,'Accelerometer');

s.DurationInSeconds = time;

s.Rate = rate;

s.Channels(1).Sensitivity = 1;
% Note:  This is not the exact sensitivity of the calibrated sensor, but it
% is close.  Was something like 1.092 on the analysis sheet, but I cannot
% find it.

s.Channels(1)

[data,time] = s.startForeground;

sampleLength = length(data);

frequencyDomainData = fft(data);

twoSidedSpectrum = abs(frequencyDomainData/sampleLength);

singleSidedSpectrum = twoSidedSpectrum(1:sampleLength/2+1);

singleSidedSpectrum(2:end-1) = 2*singleSidedSpectrum(2:end-1);

frequencyAxis = rate*(0:(sampleLength/2))/sampleLength;

yAmplitude = max(singleSidedSpectrum(200:end));

figure;
title('Frequency Domain Vibration Analysis');
subplot(2,1,1);
plot(frequencyAxis,singleSidedSpectrum);
xlabel('Frequency (Hz)');
ylabel('Amplitude (G)');

subplot(2,1,2);
plot(frequencyAxis,singleSidedSpectrum);
xlim([0 1000]);
ylim([0 yAmplitude]);
xlabel('Frequency (Hz)');
ylabel('Amplitude (G)');

release(s);

end