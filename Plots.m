% Nimish Dhawan
% March 2, 2026
% Plots

% Plotting GPS data;

dir = [cd,'\Sim Results\'];

% Position 
figure()
subplot(3,1,1)
plot(tEpoch, xHist(1,:),'ko-');
grid on; hold on;
plot(tEpoch, xp, 'r');
ylabel('x [km]',Interpreter='tex'); legend('Estimate', 'Ground Truth')
theme('light'); fontname('Times New Roman'); fontsize(10,'points');

subplot(3,1,2)
plot(tEpoch, xHist(2,:),'ko-');
grid on; hold on;
plot(tEpoch, yp, 'r');
ylabel('y [km]',Interpreter='tex');
theme('light'); fontname('Times New Roman'); fontsize(10,'points');

subplot(3,1,3)
plot(tEpoch, xHist(3,:), 'ko-');
grid on; hold on;
plot(tEpoch, zp, 'r');
ylabel('z [km]',Interpreter='tex'); xlabel('Time [s]');
theme('light'); fontname('Times New Roman'); fontsize(10,'points');
% exportgraphics(gcf, [dir, 'XYZcomparison.pdf'])

% Filter performance
figure()
subplot(3,1,1)
plot(tEpoch, xHist(1,:)-xp,'ko-');
grid on; hold on;
ylabel('\delta x [km]'); 
theme('light'); fontname('Times New Roman'); fontsize(10,'points');

subplot(3,1,2)
plot(tEpoch, xHist(2,:)-yp,'ko-');
grid on; hold on;
ylabel('\delta y [km]');
theme('light'); fontname('Times New Roman'); fontsize(10,'points');

subplot(3,1,3)
plot(tEpoch, xHist(3,:)-zp, 'ko-');
grid on; hold on;
ylabel('\delta z [km]'); xlabel('Time [s]');
theme('light'); fontname('Times New Roman'); fontsize(10,'points');
% exportgraphics(gcf, [dir, 'XYZerr.pdf'])

% Plotting covariance
figure()
grid on; hold on; box on;
for l = 1:7
    for m = 1:7
        plot(tEpoch, squeeze(PHist(l,m,:)));
    end
end
ylabel('Covariance'); xlabel('Time [s]');
theme('light'); fontname('Times New Roman'); fontsize(10,'points');
% exportgraphics(gcf, [dir, 'Covariance.pdf'])

% Plotting altitude
figure()
plot(tEpoch, altEst_km, 'ko-')
hold on; grid on; box on;
plot(tEpoch, altTruth_km, 'r-')
xlabel('Time [s]'); ylabel('Altitude [km]');
legend('Estimate', 'Ground Truth')
theme('light'); fontname('Times New Roman'); fontsize(10,'points');

% exportgraphics(gcf, [dir, 'Alt.pdf'])
