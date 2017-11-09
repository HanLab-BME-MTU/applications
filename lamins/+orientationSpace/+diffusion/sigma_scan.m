sdv = [];
A = 0.001;
B = 1;
syms sx
assume(sx,'real')

ssigma_scan_range = 0.1:0.01:1;
sdv = zeros(1,length(ssigma_scan_range));

for ssigma = ssigma_scan_range
f = A.*exp(-(sx-1).^2/2/ssigma.^2) + B.*exp(-(sx+1).^2/2/ssigma.^2);
extremum = vpasolve(diff(f,sx,1),sx);
sdv(end+1) = vpa(subs(diff(f,sx,2),sx,extremum));
end
figure; plot(ssigma_scan_range,sdv)
grid on;
xlabel('Sigma')
ylabel('Second Derivative');
title(sprintf('A = %02d, B = %d',A,B));