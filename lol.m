clear all
clc; clf;
fuckme

subplot(411);
hold on;
grid on;
title('System states');
xlabel('Time (s)');
ylabel('State (rad, ^r^a^d/_s)');
xlim([0 time(end)]);
plot(time, [ret.ssX.xPos; ret.ssX.xVel; r]');
legend('pos','vel', 'position reference');

subplot(412);
hold on;
grid on;
title('Tracking errors');
xlabel('Time (s)');
ylabel('Error (rad)');
xlim([0 time(end)]);
plot(time, [abs(ret.ssX.xPos - r); abs(ret.ssX.xPos - ret.ssL1.l1sXhat.xPos)]');
legend('position err from reference input', ...
       'position err from desired system response');

subplot(413);
hold on;
grid on;
title('Control input');
xlabel('Time (s)');
ylabel('Input torque (N-m)');
xlim([0 time(end)]);
plot(time, [ret.ssL1.l1sU; 8*ones(size(time)); -8*ones(size(time));]');
legend('desired input', 'maximum torque bound', 'minimum torque bound');

subplot(414);
hold on;
grid on;
title('Parameter estimates');
xlabel('Time (s)');
ylabel('Coefficient value');
xlim([0 time(end)]);
plot(time, [ret.ssL1.l1sWqsHat.wqsOmega; ...
            ret.ssL1.l1sWqsHat.wqsSigma; ...
            ret.ssL1.l1sWqsHat.wqsTheta.xPos; ...
            ret.ssL1.l1sWqsHat.wqsTheta.xVel; ...
           ]');
legend('omega','sigma','thetaPos','thetaVel');
