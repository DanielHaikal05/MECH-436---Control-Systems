clear all; close all;

% Experimental Bode plot
samples = 50;
W=logspace(-6, 0, samples);
mag = zeros(1,samples); ang = zeros(1,samples);
mag_pade = zeros(1,samples); ang_pade = zeros(1,samples);
[num_pade, den_pade] = pade(9.3, 1);

%Setup for CPU-parallellized 'parsim' function
model1 = 'system_id_simulink';
model2 = 'system_id_simulink';
load_system(model1)
load_system(model2)
in1(1:samples) = Simulink.SimulationInput(model1);
in2(1:samples) = Simulink.SimulationInput(model2);

for i = 1:samples
    w = W(i); 
    dt = 2*pi/(w*40);
    Tstop = max(1000, 4*2*pi/w);

    in1(i) = in1(i).setModelParameter( ...
        'StopTime', num2str(Tstop), ...
        'MaxStep', num2str(dt));
    in1(i) = in1(i).setVariable('w', w);
    in1(i) = in1(i).setVariable('use_pade', false);

    in2(i) = in2(i).setModelParameter( ...
        'StopTime', num2str(Tstop), ...
        'MaxStep', num2str(dt));
    in2(i) = in2(i).setVariable('w', w);
    in2(i) = in2(i).setVariable('use_pade', true);
    in2(i) = in2(i).setVariable('num_pade', num_pade);
    in2(i) = in2(i).setVariable('den_pade', den_pade);
end

out1 = parsim(in1, 'ShowProgress', 'on');
out2 = parsim(in2, 'ShowProgress', 'on');


%Compute magnitude and phase angle for each frequency for both true TF, and
%pade-approximated TF
for i=(1:samples)
    w=W(i);
    % Exact time delay
    y=out1(i).y.data; x=out1(i).x.data; t=out1(i).tout;
    i_ss = floor(0.7*numel(t)); i_max = numel(t);
    y_ss = y(i_ss : i_max);
    t_ss = t(i_ss : i_max);
    
    % sin(omega.t+phi) = cos(phi)sin(omega.t) + sin(phi)cos(omega.t)
    A = [sin(w*t_ss) cos(w*t_ss) ones(size(t_ss))];
    c = A \ y_ss;
    mag(i) = norm(c(1:2));
    ang(i) = atan2(c(2), c(1));

    % Pade approximation
    y=out2(i).y.data; x=out2(i).x.data; t=out2(i).tout;
    i_ss = floor(0.7*numel(t)); i_max = numel(t);
    y_ss = y(i_ss : i_max);
    t_ss = t(i_ss : i_max);
        
    % sin(omega.t+phi) = cos(phi)sin(omega.t) + sin(phi)cos(omega.t)
    A = [sin(w*t_ss) cos(w*t_ss) ones(size(t_ss))];
    c = A \ y_ss;
    mag(i) = norm(c(1:2));
    ang(i) = atan2(c(2), c(1));
    
    fprintf("Done %d/%d\n",i,samples)
end


% Transfer function obtained by system id
K_id = 2.08e-4;
z_id = [3.5e-5 1.38e-3];
p_id = [5.94e-5 3.67e-4 6.48e-3 4.33e-2];
G_id = K_id*tf(poly(-z_id),poly(-p_id));
[mag_id, phase_id, wout_id] = bode(G_id,W);
mag_id = squeeze(mag_id); phase_id = squeeze(phase_id);
Td_id = 10.36;

% Actual Bode plot
true_num = [1 1.3893e-3 5.3055e-8];
true_den = [1 3.87704e-2 2.551212584e-4 8.7843039904e-8 5.619890472e-12];
K = 1.75e-4;
G_true=K*tf(true_num, true_den, 'InputDelay', 9.3);
[true_mag, true_phase, true_w] = bode(G_true, W);
true_mag = squeeze(true_mag); true_phase = squeeze(true_phase);


% Asymptotic plot calculation
slopes = zeros(1, numel(z_id)+numel(p_id)+1);
idx_z = 1; idx_p = 1;
for i=(2:numel(slopes))
    if idx_p > numel(p_id) || (idx_z<=numel(z_id) && z_id(idx_z) < p_id(idx_p))
        slopes(i)=slopes(i-1)+20;
        idx_z = idx_z + 1;
    elseif idx_z > numel(z_id) || (idx_p<=numel(p_id) && p_id(idx_p) < z_id(idx_z))
        slopes(i)=slopes(i-1)-20;
        idx_p = idx_p + 1;
    end
end
x_approx = sort([1e-6 z_id p_id 1e-0]);
y_approx = zeros(size(x_approx));
y_approx(1) = 4.4;
for i=(2:numel(y_approx))
    dx = log10(x_approx(i)) - log10(x_approx(i-1));
    y_approx(i) = y_approx(i-1) + dx*slopes(i-1);
end


% Plotting
ang=unwrap(ang);
ang_pade=unwrap(ang_pade);

figure
semilogx(W,20*log10(mag)); hold on
semilogx(W,20*log10(mag_id)); hold on
semilogx(W,20*log10(mag_pade)); hold on
semilogx(W, 20*log10(true_mag)); hold on
semilogx(x_approx, y_approx, '--'); hold on
hold off
title('Magnitude')
legend('G_{exp}','G_{id}', 'G_{exp,pade}', 'G_{true}','Asymptotes')

figure
semilogx(W, ang*180/pi); hold on
semilogx(wout_id, phase_id); hold on
semilogx(W, ang_pade*180/pi); hold on
semilogx(W, true_phase); hold off
title('Phase')
legend('G_{exp}','G_{id}', 'G_{exp,pade}', 'G_{true}')