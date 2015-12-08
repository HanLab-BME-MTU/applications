clc;

clear;

devices=daq.getDevices;

session=daq.createSession('ni');a

session.Rate = 10000;

session.DurationInSeconds = 2;

addAnalogInputChannel(session,'Dev2',0,'Voltage');
addAnalogOutputChannel(session,'Dev2',0,'Voltage');

outputSignalCh0 = 2+2.*square(linspace(0,pi*100,session.Rate)');

outputSignalCh0(end)=0;

queueOutputData(session, outputSignalCh0);

[captured_data, time]=session.startForeground();

figure; 

subplot(2,1,1);
plot(time, abs(captured_data)./max(abs(captured_data)));

hold on;

plot(time, outputSignalCh0./max(outputSignalCh0));

hold off;

test=abs(captured_data)./max(abs(captured_data))

i=test<0.80;
test(i)=[];

subplot(2,1,2);
hist(test,200);
mean(test)
std(test)
std(test)/mean(test)