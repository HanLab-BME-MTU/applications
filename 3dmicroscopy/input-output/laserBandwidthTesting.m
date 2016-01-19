% Data Acquisition for Bandwidth Testing
% 2014-11-24, Kevin Dean

clc;
clear;

%Identify devices.
devices=daq.getDevices();
s = daq.createSession('ni');
s.Rate=1E6;
addAnalogInputChannel(s,'Dev1',0,'Voltage');
addAnalogOutputChannel(s,'Dev1',0,'Voltage');

loopNumber=0;

tic

%for i=1:1:50
    
    frequency = 100;
    
    loopNumber=loopNumber+1;
    
    %The output should always have a low voltage of -0.5, but the maximum
    %value scales linearly with the maximum voltage.
    outputSignal=1.5*square(linspace(0,pi*2*frequency,s.Rate))+1;
    outputSignal=outputSignal';
    queueOutputData(s,outputSignal);
    [data,time] = s.startForeground();
    
    queueOutputData(s,-0.5);
    s.startForeground();
    plot(time,data)
    
    %%
    freqData=fft(data);
    freqDataAmp = abs(freqData);
    
    signalAmplitude = freqDataAmp(frequency+1);
    signalAmplitude=signalAmplitude*2;
    signalAmplitude=signalAmplitude/1E5;
    
    freqDataPhase = unwrap(angle(freqData));
    signalPhase = freqDataPhase(frequency+1);
    
    FinalData(loopNumber,1)=frequency;
    FinalData(loopNumber,2)=signalPhase;
    FinalData(loopNumber,3)=signalAmplitude;
    
end

toc

subplot(2,1,1);
plot(FinalData(:,1),FinalData(:,3),'x-');

subplot(2,1,2);
plot(FinalData(:,1),FinalData(:,2),'x-');


%%


subplot(2,1,1);
plot(water(:,1),water(:,3),'kx-');
hold on;
plot(air(:,1),air(:,3),'o--r')
hold off

subplot(2,1,2);
plot(water(:,1),water(:,2),'x-');
hold on;
plot(air(:,1),air(:,2),'o--r')
hold off
