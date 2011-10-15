

%dt = 0.1;
%t = 0:dt:50;
% S_init = [1 0];
% kVect = 1;
% [ti, Y] = ode45(@(t,y) singleExpModel(t, y, kVect), [t(1) t(end)], S_init);
% figure;
% plot(ti, Y);


% Cannot work:
% S_init = [1 0 0];
% kVect = [1 0.5];
% [ti, Y] = ode45(@(t,y) doubleExpModel(t, y, kVect), [t(1) t(end)], S_init);
% figure;
% plot(ti, Y);
% hold on;
% plot(ti, Y(:,2)/Y(end,2), 'm');
% plot(ti, Y(:,3)/Y(end,3), 'c--');

%plot(t, Y(end,2)*exp(-kVect(1)*t) + Y(end,3)*exp(-kVect(2)*t), 'r--');
% plot(t, Y(end,2)*exp(-kVect(1)*t), 'r--');
% plot(t, Y(end,2)-Y(end,2)*exp(-1.5*t), 'm--');
% plot(t, Y(end,3)*exp(-1.5*t), 'c--');




% % dt = t(2)-t(1);
% % 
% % S_init = [1 0 0 0];
% % kVect = [0.5 0.5 0.3];
% % [ti, Y] = ode45(@(t,y) pop2Model(t, y, kVect), [t(1) t(end)], S_init);
% % % figure;
% % % plot(ti, Y);
% % % hold on;
% % % plot(ti, sum(Y(:,[1 3]), 2), 'k--');
% % 
% % 
% % 
% % % output weights
% % w = Y(end, [2 4]);
% % 
% % % PDF
% % p1 = interp1(ti, Y(:,1)/sum(Y(:,1)*dt), t);
% % p2 = interp1(ti, Y(:,3)/sum(Y(:,3)*dt), t);
% % pdf = w(1)*p1 + w(2)*p2;
% % %figure;
% % %plot(t, pdf);



% figure;
% k1 = 0.9;
% k2 = 0.1;
% f = -k1*k2*(exp(-k1*t)-exp(-k2*t))/(k1-k2);
% hold on;
% plot(t, f, 'm--');
% 
% 
% sum(f*dt)



%  function v = singleExpCost(kVect, signal, t, S_init)
% kVect = abs(kVect);
% [ti, Y] = ode45(@(t,y) modelC(t, y, kVect), [2 t(end)], S_init);
% %model = Y(:,2)*Y(end,6) + Y(:,3)*Y(end,7) + Y(:,4)*Y(end,5);
% %model = model / sum(model.*[0; diff(ti)]);
% model = sum(Y(:,[5 6 7]),2);
% 
% v = signal - interp1(ti, model, t);


% function v = doubleExpCost(p, t, f_hist)
% 
% S_init = [1 0 0];
% [ti, Y] = ode45(@(t,y) doubleExpModel(t, y, p), [t(1) t(end)], S_init);
% 
% figure;
% plot(ti, Y);
% 


% function dy = singleExpModel(~, y, k)
% N = numel(y);
% dy = zeros(N,1);
% dy(1) = -k(1)*y(1);
% dy(2) = k(1)*y(1);
% 
% 
% function dy = doubleExpModel(~, y, k)
% N = numel(y);
% dy = zeros(N,1);
% dy(1) = -(k(1)+k(2))*y(1);
% dy(2) = k(1)*y(1);
% dy(3) = k(2)*y(1);





%  function v = costC(kVect, signal, t, S_init)
% kVect = abs(kVect);
% [ti, Y] = ode45(@(t,y) modelC(t, y, kVect), [2 t(end)], S_init);
% %model = Y(:,2)*Y(end,6) + Y(:,3)*Y(end,7) + Y(:,4)*Y(end,5);
% %model = model / sum(model.*[0; diff(ti)]);
% model = sum(Y(:,[5 6 7]),2);
% 
% v = signal - interp1(ti, model, t);
% 
% function dy = modelC(~, y, k)
% N = length(y);
% dy = zeros(N,1);
% dy(1) = -k(1)*y(1);
% dy(2) = -(k(2)+k(5))*y(2) + k(1)*y(1);
% dy(3) = -(k(3)+k(6))*y(3) + k(2)*y(2);
% dy(4) = -k(4)*y(4) + k(3)*y(3);
% dy(5) = k(4)*y(4);
% dy(6) = k(5)*y(2);
% dy(7) = k(6)*y(3);




% dt = [0 diff(t)];
% EDF = cumsum(data.*dt);
% EDF = EDF - max(EDF) + 1;
% 
% 
% %===================================
% % Model C
% %===================================
% 
% % data plot
% figure;
% set(gcf, 'Position',  [440 358 560 320], 'PaperPositionMode', 'auto');
% plot(t, EDF, 'k.', 'MarkerSize', 12);
% %axis([0 100 0 0.1]);
% axis([0 t(end) 0 1]);
% hold on;
% 
% 
% % initialization
% nStates = 7;
% S_init = zeros(1, nStates);
% S_init(1) = 1;
% 
% kInit = 0.1*rand(1,nStates);
% kInit(1) = 1;
% 
% %kInit = [0.8140    0.0753    0.0233    0.0220    0.3250    0.0556    0.0656];
% % kInit = [0.5104    0.3635    0.0375    0.0146    0.1503    0.0325    0.0853];
% 
% [kVect resnorm] = lsqnonlin(@costC, kInit, [], [], opts, EDF, t, S_init);
% kVect = abs(kVect)
% 
% n = length(data); % data points
% deg = length(kVect); % degrees of freedom
% BIC = n*log(resnorm/n) + deg*log(n);
% 
% 
% [ti,Y_fit] = ode45(@(t0,y0) modelC(t0, y0, kVect), [2 t(end)], S_init);
% 
% 
% modelPDF = Y_fit(:,2)*Y_fit(end,6) + Y_fit(:,3)*Y_fit(end,7) + Y_fit(:,4)*Y_fit(end,5);
% modelCDF = sum(Y_fit(:,[5 6 7]),2);
% 
% res.A = [Y_fit(end,6) Y_fit(end,7) Y_fit(end,5)]
% 
% 
% plot(ti, modelCDF, 'r');
% 
% % figure;
% % set(gcf, 'Position',  [440 358 560 320], 'PaperPositionMode', 'auto');
% % plot(t, data, 'k.', 'MarkerSize', 12);
% % axis([0 100 0 0.1]);
% % hold on;
% % plot(ti, modelPDF/3, 'r');
% 
% figure;
% plot(t, [0 diff(EDF)]);
% hold on;
% plot(ti, [0; diff(modelCDF)], 'r');
% 
% 
% %plot(ti, Y_fit);
% return
% % model = Y_fit(:,2)*Y_fit(end,6) + Y_fit(:,3)*Y_fit(end,7) + Y_fit(:,4)*Y_fit(end,5);
% % nrm = sum(model.*[0; diff(ti)]);
% % model = model / nrm;
% 
% 
% P1 = Y_fit(:,2)/nrm;
% P2 = Y_fit(:,3)/nrm;
% P3 = Y_fit(:,4)/nrm;
% 
% dt = [diff(ti); 0];
% p1 = P1/sum(P1.*dt);
% p2 = P2/sum(P2.*dt);
% p3 = P3/sum(P3.*dt);
% 
% tau1 = sum(p1.*ti.*dt);
% tau2 = sum(p2.*ti.*dt);
% tau3 = sum(p3.*ti.*dt);
% 
% [cs1 idx1] = unique(cumsum(p1.*dt));
% [cs2 idx2] = unique(cumsum(p2.*dt));
% [cs3 idx3] = unique(cumsum(p3.*dt));
% 
% res.p25 = [interp1(cs1, ti(idx1), 0.25) interp1(cs2, ti(idx2), 0.25) interp1(cs3, ti(idx3), 0.25)];
% res.p75 = [interp1(cs1, ti(idx1), 0.75) interp1(cs2, ti(idx2), 0.75) interp1(cs3, ti(idx3), 0.75)];
% 
% res.A = [Y_fit(end,6) Y_fit(end,7) Y_fit(end,5)];
% res.tau = [tau1 tau2 tau3];
% 
% fprintf('SE: %s, BIC: %.2f\n', resnorm, BIC);
% fprintf('Amplitudes: A1 = %2.2f, A2 = %2.2f, A3 = %2.2f\n', res.A(1), res.A(2), res.A(3));
% fprintf('Means: tau1 = %2.2f, tau2 = %2.2f, tau3 = %2.2f\n', tau1, tau2, tau3);
% 
% 
% model = interp1(ti, model, t);
% Y_fit = interp1(ti, Y_fit, t);
% 
% plot(t, model, '--', 'Color', [1 0 0], 'LineWidth', 2);
% %plot(t, Y_fit(:,[2 3 4])/nrm, '-', 'Color', [0.6 0 0], 'LineWidth', 1);
% 
% plot(t, res.A(1)*Y_fit(:,2)/nrm, '-', 'Color', [0.6 0 0], 'LineWidth', 1);
% plot(t, res.A(2)*Y_fit(:,3)/nrm, '-', 'Color', [0.6 0 0], 'LineWidth', 1);
% plot(t, res.A(3)*Y_fit(:,4)/nrm, '-', 'Color', [0.6 0 0], 'LineWidth', 1);
% 
% set(gca, 'FontName', 'Helvetica', 'FontSize', 14, 'LineWidth', 1.5);
% xlabel('Time (s)', 'FontName', 'Helvetica', 'FontSize', 14);
% ylabel('Relative frequency', 'FontName', 'Helvetica', 'FontSize', 14);
% legend('Lifetime distribution', 'Model', 'State lifetimes', 'Location', 'NorthEast');
% %print('-depsc2', '-r300', 'modelExample.eps');
% 
% % plot state outputs
% figure;
% set(gcf, 'Position',  [440 358 560 320], 'PaperPositionMode', 'auto');
% plot(t, Y_fit(:,[5 6 7]), '-', 'Color', [0.6 0 0], 'LineWidth', 1);
% %hold on;
% %plot(t, Y_fit(:,[2 3 4]), '-', 'Color', [0 1 0], 'LineWidth', 1);
% %plot(t, Y_fit(:,1), 'k-', 'LineWidth', 1);
% set(gca, 'FontName', 'Helvetica', 'FontSize', 14, 'LineWidth', 1.5);
% xlabel('Time (s)', 'FontName', 'Helvetica', 'FontSize', 14);
% ylabel('Relative frequency', 'FontName', 'Helvetica', 'FontSize', 14);
% axis([0 t(end) 0 0.7]);
% 
% 
% 
% 
% function v = costC(kVect, signal, t, S_init)
% kVect = abs(kVect);
% [ti, Y] = ode45(@(t,y) modelC(t, y, kVect), [2 t(end)], S_init);
% %model = Y(:,2)*Y(end,6) + Y(:,3)*Y(end,7) + Y(:,4)*Y(end,5);
% %model = model / sum(model.*[0; diff(ti)]);
% model = sum(Y(:,[5 6 7]),2);
% 
% v = signal - interp1(ti, model, t);
% 
% function dy = modelC(~, y, k)
% N = length(y);
% dy = zeros(N,1);
% dy(1) = -k(1)*y(1);
% dy(2) = -(k(2)+k(5))*y(2) + k(1)*y(1);
% dy(3) = -(k(3)+k(6))*y(3) + k(2)*y(2);
% dy(4) = -k(4)*y(4) + k(3)*y(3);
% dy(5) = k(4)*y(4);
% dy(6) = k(5)*y(2);
% dy(7) = k(6)*y(3);
% 
% 
