clear; clc; close all;
set(0,'DefaultFigureVisible','on');
s=tf('s');

% PART III PLANT
K=2.33e-1;
z1=6.12e-2;
z2=8.55e-2;
p1=9.71e-2;
p2=2.68e-2;
p3=4.34e-3;
p4=0.52e-3;
Td=11.7;
G_nodelay=K*(s+z1)*(s+z2)/((s+p1)*(s+p2)*(s+p3)*(s+p4));
% Exact delayed plant for verification
G_exact=zpk([-z1 -z2],[-p1 -p2 -p3 -p4],K,'InputDelay',Td);
% First-order Pade plant for Bode/margins/design checking
[numD,denD]=pade(Td,1);
G_pade=minreal(tf(numD,denD)*G_nodelay);

% TIME VECTOR
t=0:1:15000;
% SPECS AND IDEAL SECOND-ORDER SYSTEM
Mp=0.03;
Ts_required=1200;
ess_required=0.02;
zeta=-log(Mp)/sqrt(pi^2+log(Mp)^2);
wn=4/(zeta*Ts_required);
Tideal=tf(wn^2,[1 2*zeta*wn wn^2]);
[yid,tid]=step(Tideal,t);
yid=squeeze(yid);
fprintf('\nIDEAL SECOND ORDER\n');
fprintf('zeta=%.6f\n',zeta);
fprintf('wn=%.6e rad/s\n',wn);

% NEW FREQUENCY-RESPONSE LEAD-LAG CONTROLLER
numLL=[1.73797606e-04 5.57358725e-07 2.53973523e-10];
denLL=[1 1.14727809e-02 9.96311098e-07];
C_LL=tf(numLL,denLL);
fprintf('\nNEW FREQUENCY-RESPONSE LEAD-LAG CONTROLLER\n');
fprintf('Numerator:\n'); disp(numLL);
fprintf('Denominator:\n'); disp(denLL);

% NEW FREQUENCY-RESPONSE PID CONTROLLER
numPID=[3.93864393e-04 2.31416100e-06 1.17709441e-09];
denPID=[1 4.50000000e-02 0];
C_PID=tf(numPID,denPID);
% Equivalent PID block values
Kp=5.084451890329e-05;
Ki=2.615765353720e-08;
Kd=7.622663874609e-03;
N=4.500000000000e-02;
fprintf('\nNEW FREQUENCY-RESPONSE PID CONTROLLER\n');
fprintf('Transfer function numerator:\n'); disp(numPID);
fprintf('Transfer function denominator:\n'); disp(denPID);
fprintf('\nEquivalent PID block parameters:\n');
fprintf('Kp=%.12e\n',Kp);
fprintf('Ki=%.12e\n',Ki);
fprintf('Kd=%.12e\n',Kd);
fprintf('N=%.12e\n',N);

% The PID controller is also defined from frequency-response design.
% The equivalent Kp, Ki, Kd, and filter N are extracted for Simulink use.

% CLOSED-LOOP SYSTEMS ON ORIGINAL DELAYED PLANT
T_LL=feedback(C_LL*G_exact,1);
T_PID=feedback(C_PID*G_exact,1);
[yLL,tLL]=step(T_LL,t);
[yPID,tPID]=step(T_PID,t);
yLL=squeeze(yLL);
yPID=squeeze(yPID);
info_LL=stepinfo(yLL,tLL,dcgain(T_LL));
info_PID=stepinfo(yPID,tPID,dcgain(T_PID));
ess_LL=abs(1-dcgain(T_LL));
ess_PID=abs(1-dcgain(T_PID));
fprintf('\nLEAD-LAG PERFORMANCE\n');
fprintf('Ts=%.2f s\n',info_LL.SettlingTime);
fprintf('OS=%.2f %%\n',info_LL.Overshoot);
fprintf('ess=%.5f\n',ess_LL);
fprintf('\nPID PERFORMANCE\n');
fprintf('Ts=%.2f s\n',info_PID.SettlingTime);
fprintf('OS=%.2f %%\n',info_PID.Overshoot);
fprintf('ess=%.5f\n',ess_PID);

% Both controllers are tested on the exact delayed plant to evaluate real performance.
% We compare settling time, overshoot, and steady-state error.

fprintf('\nSPEC CHECK\n');
check_specs(info_LL,ess_LL,'Lead-Lag');
check_specs(info_PID,ess_PID,'PID');

% CONTROL EFFORT USING PADE MODEL
U_LL=minreal(C_LL/(1+C_LL*G_pade));
U_PID=minreal(C_PID/(1+C_PID*G_pade));
[uLL,tuLL]=step(U_LL,t);
[uPID,tuPID]=step(U_PID,t);
uLL=squeeze(uLL);
uPID=squeeze(uPID);
fprintf('\nCONTROL EFFORT\n');
fprintf('Lead-Lag umax=%.6e\n',max(abs(uLL)));
fprintf('PID umax=%.6e\n',max(abs(uPID)));

% This compares how aggressive each controller is in terms of control input.
% Higher values mean more effort and possibly less practical implementation.

% FREQUENCY MARGINS USING GPADE
[GMll,PMll,Wcgll,Wcpll]=margin(C_LL*G_pade);
[GMpid,PMpid,Wcgpid,Wcppid]=margin(C_PID*G_pade);
fprintf('\nFREQUENCY MARGINS\n');
fprintf('Lead-Lag GM=%.2f dB | PM=%.2f deg\n',20*log10(GMll),PMll);
fprintf('PID GM=%.2f dB | PM=%.2f deg\n',20*log10(GMpid),PMpid);
DMll=deg2rad(PMll)/Wcpll;
DMpid=deg2rad(PMpid)/Wcppid;
fprintf('Lead-Lag delay margin approx=%.2f s\n',DMll);
fprintf('PID delay margin approx=%.2f s\n',DMpid);

% Gain margin, phase margin, and delay margin give insight about robustness.
% They show how much uncertainty or delay the system can tolerate.

% DELAY ROBUSTNESS
delay_factors=[1 1.2 1.4 1.6];
fprintf('\nDELAY ROBUSTNESS\n');
YLL_delay=cell(size(delay_factors));
YPID_delay=cell(size(delay_factors));
Tdelay_store=cell(size(delay_factors));
for k = 1:length(delay_factors)

    Td_new=Td*delay_factors(k);
    G_var=zpk([-z1 -z2],[-p1 -p2 -p3 -p4],K,'InputDelay',Td_new);
    T_LL_var=feedback(C_LL*G_var,1);
    T_PID_var=feedback(C_PID*G_var,1);
    [yLL_var,tvar]=step(T_LL_var,t);
    [yPID_var,~]=step(T_PID_var,t);
    yLL_var=squeeze(yLL_var);
    yPID_var=squeeze(yPID_var);
    YLL_delay{k}=yLL_var;
    YPID_delay{k}=yPID_var;
    Tdelay_store{k}=tvar;
    info_LL_var=stepinfo(yLL_var,tvar,dcgain(T_LL_var));
    info_PID_var=stepinfo(yPID_var,tvar,dcgain(T_PID_var));
    fprintf('\nDelay factor=%.1f\n',delay_factors(k));
    fprintf('Lead-Lag Ts=%.2f | OS=%.2f\n', ...
        info_LL_var.SettlingTime, info_LL_var.Overshoot);
    fprintf('PID Ts=%.2f | OS=%.2f\n', ...
        info_PID_var.SettlingTime, info_PID_var.Overshoot);
end

% This part studies how performance changes when the delay increases slightly.
% It helps evaluate robustness of both controllers to modeling errors.

% EXTREME DELAY ROBUSTNESS: 10Td, 20Td, 30Td
delay_factors_extreme=[10 20 30];
t_extreme=0:1:60000;
fprintf('\nEXTREME DELAY ROBUSTNESS: 10Td, 20Td, 30Td\n');
YLL_extreme=cell(size(delay_factors_extreme));
YPID_extreme=cell(size(delay_factors_extreme));
Textreme_store=cell(size(delay_factors_extreme));
Ts_LL_extreme=nan(size(delay_factors_extreme));
OS_LL_extreme=nan(size(delay_factors_extreme));
ess_LL_extreme=nan(size(delay_factors_extreme));
Ts_PID_extreme=nan(size(delay_factors_extreme));
OS_PID_extreme=nan(size(delay_factors_extreme));
ess_PID_extreme=nan(size(delay_factors_extreme));
for k=1:length(delay_factors_extreme)
    Td_new = Td * delay_factors_extreme(k);
    G_var=zpk([-z1 -z2],[-p1 -p2 -p3 -p4],K,'InputDelay',Td_new);
    T_LL_var=feedback(C_LL*G_var,1);
    T_PID_var=feedback(C_PID*G_var,1);
    [yLL_var,tLL_var]=step(T_LL_var,t_extreme);
    [yPID_var,tPID_var]=step(T_PID_var,t_extreme);
    yLL_var=squeeze(yLL_var);
    yPID_var=squeeze(yPID_var);
    YLL_extreme{k}=yLL_var;
    YPID_extreme{k}=yPID_var;
    Textreme_store{k}=tLL_var;
    % Final values
    yfinal_LL=dcgain(T_LL_var);
    yfinal_PID=dcgain(T_PID_var);
    % Step info
    info_LL_var=stepinfo(yLL_var,tLL_var,yfinal_LL);
    info_PID_var=stepinfo(yPID_var,tPID_var,yfinal_PID);
    % Store metrics
    Ts_LL_extreme(k)=info_LL_var.SettlingTime;
    OS_LL_extreme(k)=info_LL_var.Overshoot;
    ess_LL_extreme(k)=abs(1-yfinal_LL);
    Ts_PID_extreme(k)=info_PID_var.SettlingTime;
    OS_PID_extreme(k)=info_PID_var.Overshoot;
    ess_PID_extreme(k)=abs(1 - yfinal_PID);
    fprintf('\nDelay=%dTd\n',delay_factors_extreme(k));
    fprintf('Lead-Lag Ts=%.2f s | OS=%.2f %% | ess=%.5f\n', ...
        Ts_LL_extreme(k), OS_LL_extreme(k), ess_LL_extreme(k));

    fprintf('PID Ts=%.2f s | OS=%.2f %% | ess=%.5f\n', ...
        Ts_PID_extreme(k), OS_PID_extreme(k), ess_PID_extreme(k));
end

% Extreme delay cases (10Td, 20Td, 30Td) show the limits of stability.
% If responses become very large or unstable, it means the controller is not robust.

% Summary table for extreme delay cases
ExtremeDelayTable=table( ...
    delay_factors_extreme', ...
    Ts_LL_extreme', OS_LL_extreme', ess_LL_extreme', ...
    Ts_PID_extreme', OS_PID_extreme', ess_PID_extreme', ...
    'VariableNames', {'DelayFactor','Ts_LL','OS_LL','ess_LL','Ts_PID','OS_PID','ess_PID'});
disp(' ');
disp('EXTREME DELAY SUMMARY TABLE');
disp(ExtremeDelayTable);

% MATLAB PLOTS
figure;
plot(tLL,yLL,'b','LineWidth',1.5);
grid on;
xlabel('Time [s]');
ylabel('Output');
title('New Frequency-Response Lead-Lag Step Response');
figure;
plot(tPID,yPID,'r','LineWidth',1.5);
grid on;
xlabel('Time [s]');
ylabel('Output');
title('New Frequency-Response PID Step Response');
figure;
plot(tLL,yLL,'b','LineWidth',1.5); hold on;
plot(tPID,yPID,'r','LineWidth',1.5);
plot(tid,yid,'k--','LineWidth',1.5);
legend('Lead-Lag','PID','Ideal 2nd Order','Location','best');
grid on;
xlabel('Time [s]');
ylabel('Output');
title('New Frequency-Response Closed-Loop Response Comparison');
figure;
plot(tuLL,uLL,'b','LineWidth',1.5); hold on;
plot(tuPID,uPID,'r','LineWidth',1.5);
legend('Lead-Lag','PID','Location','best');
grid on;
xlabel('Time [s]');
ylabel('u(t)');
title('New Frequency-Response Control Effort');
figure;
margin(C_LL*G_pade);
grid on;
title('New Frequency-Response Lead-Lag Open-Loop Margin');
figure;
margin(C_PID*G_pade);
grid on;
title('New Frequency-Response PID Open-Loop Margin');
for k = 1:length(delay_factors)
    figure;
    plot(Tdelay_store{k},YLL_delay{k},'b','LineWidth',1.5);
    grid on;
    xlabel('Time [s]');
    ylabel('Output');
    title(['New Frequency-Response Lead-Lag Delay x',num2str(delay_factors(k))]);
    figure;
    plot(Tdelay_store{k},YPID_delay{k},'r','LineWidth',1.5);
    grid on;
    xlabel('Time [s]');
    ylabel('Output');
    title(['New Frequency-Response PID Delay x',num2str(delay_factors(k))]);
end

% Extreme delay plots: individual figures
for k = 1:length(delay_factors_extreme)
    figure;
    plot(Textreme_store{k},YLL_extreme{k},'b','LineWidth',1.5);
    grid on;
    xlabel('Time [s]');
    ylabel('Output');
    title(['Frequency-Response Lead-Lag Extreme Delay ',num2str(delay_factors_extreme(k)),'T_d']);
    figure;
    plot(Textreme_store{k},YPID_extreme{k},'r','LineWidth',1.5);
    grid on;
    xlabel('Time [s]');
    ylabel('Output');
    title(['Frequency-Response PID Extreme Delay ',num2str(delay_factors_extreme(k)),'T_d']);
end

% Extreme delay plots: combined figures
figure;
plot(Textreme_store{1},YLL_extreme{1},'LineWidth',1.5); hold on;
plot(Textreme_store{2},YLL_extreme{2},'LineWidth',1.5);
plot(Textreme_store{3},YLL_extreme{3},'LineWidth',1.5);
grid on;
xlabel('Time [s]');
ylabel('Output');
title('Frequency-Response Lead-Lag: Extreme Delay Cases');
legend('10T_d','20T_d','30T_d','Location','best');
figure;
plot(Textreme_store{1},YPID_extreme{1},'LineWidth',1.5); hold on;
plot(Textreme_store{2},YPID_extreme{2},'LineWidth',1.5);
plot(Textreme_store{3},YPID_extreme{3},'LineWidth',1.5);
grid on;
xlabel('Time [s]');
ylabel('Output');
title('Frequency-Response PID: Extreme Delay Cases');
legend('10T_d','20T_d','30T_d','Location','best');

% AUTO-SAVE FIGURES
outFolder=fullfile(pwd,'figuresfreq_new');
if ~exist(outFolder,'dir')
    mkdir(outFolder);
end
figHandles=findall(0,'Type','figure');
[~,idx]=sort([figHandles.Number]);
figHandles=figHandles(idx);
figureNames={
    'newfreq_base_leadlag'
    'newfreq_base_pid'
    'newfreq_closed_loop_comparison'
    'newfreq_control_effort'
    'newfreq_leadlag_open_loop_margin'
    'newfreq_pid_open_loop_margin'
    'newfreq_leadlag_delay_x1p0'
    'newfreq_pid_delay_x1p0'
    'newfreq_leadlag_delay_x1p2'
    'newfreq_pid_delay_x1p2'
    'newfreq_leadlag_delay_x1p4'
    'newfreq_pid_delay_x1p4'
    'newfreq_leadlag_delay_x1p6'
    'newfreq_pid_delay_x1p6'
};
for i = 1:min(length(figHandles),length(figureNames))
    fig=figHandles(i);
    filename=fullfile(outFolder,[figureNames{i},'.png']);
    figure(fig);
    set(fig,'Color','w');
    exportgraphics(fig,filename,'Resolution',300);
end
fprintf('\nAll MATLAB figures saved in:\n%s\n',outFolder);
if ispc
    winopen(outFolder);
end

% SIMULINK VALUES
fprintf('\nSIMULINK VALUES\n');
[numG,denG]=tfdata(G_nodelay,'v');
fprintf('\nPlant numerator:\n');
disp(numG);
fprintf('Plant denominator:\n');
disp(denG);
fprintf('\nTransport delay=%.2f s\n',Td);
fprintf('\nLead-Lag numerator:\n');
disp(numLL);
fprintf('Lead-Lag denominator:\n');
disp(denLL);
fprintf('\nPID transfer function numerator:\n');
disp(numPID);
fprintf('PID transfer function denominator:\n');
disp(denPID);
fprintf('\nPID block parameters:\n');
fprintf('Kp=%.12e\n',Kp);
fprintf('Ki=%.12e\n',Ki);
fprintf('Kd=%.12e\n',Kd);
fprintf('N=%.12e\n',N);

% LOCAL FUNCTION
function check_specs(info,ess,name)
fprintf('\n%s:\n',name);
fprintf('Ts=%.2f s\n',info.SettlingTime);
fprintf('OS=%.2f %%\n',info.Overshoot);
fprintf('ess=%.5f\n',ess);
if info.SettlingTime<=1200
    disp('Ts PASS');
else
    disp('Ts FAIL');
end
if info.Overshoot <= 3
    disp('OS PASS');
else
    disp('OS FAIL');
end
if ess<=0.02
    disp('ess PASS');
else
    disp('ess FAIL');
end
end

Td
numG
denG
numLL
denLL
numPID
denPID
Kp
Ki
Kd
N