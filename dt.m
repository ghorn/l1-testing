clear all; clc; clf;

solutions = dlmread("butts");
times = dlmread("times");
control = dlmread("ctrl")';
ref = dlmread("reference")';

figure(1);

h11 = subplot(311);
plot(times, solutions(:, 1:3), '-',
     times, ref(:, 1:2), '.');
hold on;
grid on;
xlim([0 times(end)]);
ylim([-1 1]);
title('System states');
xlabel('Time (s)');
ylabel('State');
legend('x_0', 'x_1', 'x_2',
       'ref_0', 'ref_1');

h12 = subplot(312);
plot(times, control(:, 4:5));
hold on;
grid on;
xlim([0 times(end)]);
title('Control inputs');
xlabel('Time (s)');
ylabel('Input');
legend('u_0', 'u_1');

h13 = subplot(313);
plot(times, control(:, 11:13), '--');
hold on;
grid on;
xlim([0 times(end)]);
title('Estimated parameters');
xlabel('Time (s)');
ylabel('Parameter');
legend('\sigma_m_0', '\sigma_m_1', '\sigma_u_m');

# linkaxes([h11, h12, h13], "x");
