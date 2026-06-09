clc; close all;

outFolder=fullfile(pwd, 'part4_extreme_delay_plots');
if ~exist(outFolder,'dir')
    mkdir(outFolder);
end
fprintf('\nSaving Part 4 extreme delay figures to:\n%s\n\n', outFolder);
set(0,'DefaultLineLineWidth',2)
set(0,'DefaultAxesFontSize',14)
ref=1;

% CHECK REQUIRED VARIABLES
requiredFields={ ...
    'y_ll_d10','y_pid_d10', ...
    'y_ll_d20','y_pid_d20', ...
    'y_ll_d30','y_pid_d30'};

for i=1:length(requiredFields)
    try
        temp=out.(requiredFields{i});
    catch
        error('Missing out.%s. Check the To Workspace block name in Simulink.', requiredFields{i});
    end
end

% DEFINE CASES
cases={
    'Lead-Lag 10Td',out.y_ll_d10;
    'PID 10Td',out.y_pid_d10;
    'Lead-Lag 20Td',out.y_ll_d20;
    'PID 20Td',out.y_pid_d20;
    'Lead-Lag 30Td',out.y_ll_d30;
    'PID 30Td',out.y_pid_d30;
};
n=size(cases,1);
ControllerCase=strings(n,1);
Ts=nan(n,1);
OS=nan(n,1);
ess=nan(n,1);
yFinal_est=nan(n,1);
yMax=nan(n,1);

fprintf('\nPART 4 SIMULINK EXTREME DELAY PERFORMANCE\n');

% COMPUTE Ts, OVERSHOOT, ess
for i=1:n
    ControllerCase(i)=cases{i,1};
    data=cases{i,2};
    [t,y]=getTimeAndSignal(data);
    % Estimate final value from last 5% of the response
    Nsamples=length(y);
    lastN=max(10,round(0.05*Nsamples));
    yFinal_est(i)=mean(y(end-lastN+1:end));
    % Steady-state error
    ess(i)=abs(ref-yFinal_est(i));
    % Maximum output
    yMax(i)=max(y);
    % Step response information
    try
        info=stepinfo(y,t,yFinal_est(i),'SettlingTimeThreshold',0.02);
        Ts(i)=info.SettlingTime;
        OS(i)=info.Overshoot;
    catch
        Ts(i)=NaN;
        OS(i)=NaN;
    end

    fprintf('\n%s\n',ControllerCase(i));
    fprintf('Estimated final value=%.5f\n',yFinal_est(i));
    fprintf('ess=%.5f\n',ess(i));
    fprintf('Ts=%.2f s\n',Ts(i));
    fprintf('OS=%.2f %%\n',OS(i));

end

% SUMMARY TABLE
ResultsTable=table(ControllerCase, Ts, OS, ess, yFinal_est, yMax);

disp(' ');
disp('===== PART 4 EXTREME DELAY SUMMARY TABLE =====');
disp(ResultsTable);
writetable(ResultsTable, fullfile(outFolder,'part4_extreme_delay_metrics.csv'));

% INDIVIDUAL FIGURES
saveplot(out.y_ll_d10,'Part 4 Lead-Lag Response (10Td)',...
    'Output','b',fullfile(outFolder,'part4_leadlag_10Td.png'))
saveplot(out.y_pid_d10,'Part 4 PID Response (10Td)',...
    'Output','r',fullfile(outFolder,'part4_pid_10Td.png'))
saveplot(out.y_ll_d20,'Part 4 Lead-Lag Response (20Td)',...
    'Output','b',fullfile(outFolder,'part4_leadlag_20Td.png'))
saveplot(out.y_pid_d20,'Part 4 PID Response (20Td)',...
    'Output','r',fullfile(outFolder,'part4_pid_20Td.png'))
saveplot(out.y_ll_d30,'Part 4 Lead-Lag Response (30Td)',...
    'Output','b',fullfile(outFolder,'part4_leadlag_30Td.png'))
saveplot(out.y_pid_d30,'Part 4 PID Response (30Td)',...
    'Output','r',fullfile(outFolder,'part4_pid_30Td.png'))

% COMBINED FIGURES
savecombined3(out.y_ll_d10, out.y_ll_d20, out.y_ll_d30, ...
    'Part 4 Lead-Lag Response: Extreme Delay Cases', ...
    'Output', ...
    {'10T_d','20T_d','30T_d'}, ...
    fullfile(outFolder,'part4_leadlag_extreme_delays_combined.png'))
savecombined3(out.y_pid_d10, out.y_pid_d20, out.y_pid_d30, ...
    'Part 4 PID Response: Extreme Delay Cases', ...
    'Output', ...
    {'10T_d','20T_d','30T_d'}, ...
    fullfile(outFolder,'part4_pid_extreme_delays_combined.png'))
fprintf('\nAll Part 4 plots and metrics saved successfully in:\n%s\n', outFolder);
if ispc
    winopen(outFolder);
end

% HELPER FUNCTION: EXTRACT TIME AND SIGNAL
function [t,y]=getTimeAndSignal(data)

    t=data.time;
    y=data.signals.values;

    if size(y,2)>1
        y=y(:,1);
    end

    t=t(:);
    y=y(:);

end

% HELPER FUNCTION: SINGLE PLOT
function saveplot(data, title_str, ylabel_str, color, filename)

    [t,y]=getTimeAndSignal(data);
    figure('Color','w')
    plot(t, y, color)
    xlabel('Time (s)')
    ylabel(ylabel_str)
    title(title_str)
    grid on
    ax=gca;
    ax.Toolbar.Visible='off';

    exportgraphics(gcf, filename, 'Resolution', 300)

end

% HELPER FUNCTION: COMBINED PLOT
function savecombined3(data1, data2, data3, title_str, ylabel_str, legend_str, filename)
    [t1,y1]=getTimeAndSignal(data1);
    [t2,y2]=getTimeAndSignal(data2);
    [t3,y3]=getTimeAndSignal(data3);
    figure('Color','w')
    plot(t1, y1, 'LineWidth', 2); hold on;
    plot(t2, y2, 'LineWidth', 2);
    plot(t3, y3, 'LineWidth', 2);
    xlabel('Time (s)')
    ylabel(ylabel_str)
    title(title_str)
    legend(legend_str, 'Location', 'best')
    grid on
    ax=gca;
    ax.Toolbar.Visible='off';
    exportgraphics(gcf, filename, 'Resolution', 300)

end