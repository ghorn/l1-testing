clear all
clc; clf;
fuckme

subplot(411)
plot(time, [ret.ssX.xPos; ret.ssX.xVel]','.-')
legend('pos','vel')

subplot(412)
plot(time, [ret.ssX.xPos - ret.ssL1.l1sXhat.xPos; ret.ssX.xVel - ret.ssL1.l1sXhat.xVel]','.-')
legend('pos err','vel err')

subplot(413)
plot(time, [ret.ssL1.l1sU]','.-')
legend('u')

subplot(414)
plot(time, [ret.ssL1.l1sWqsHat.wqsOmega]','.-')
hold on
plot(time, [ret.ssL1.l1sWqsHat.wqsSigma]','.-')
plot(time, [ret.ssL1.l1sWqsHat.wqsTheta.xPos]','.-')
plot(time, [ret.ssL1.l1sWqsHat.wqsTheta.xVel]','.-')
legend('omega','sigma','thetaPos','thetaVel')