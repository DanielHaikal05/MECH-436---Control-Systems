clear all; close all;

% Experimental Bode plot
samples = 40;
W=logspace(-5, -1, samples);
mag = zeros(1,samples); ang = zeros(1,samples);
mag_noisy = zeros(1,samples); ang_noisy = zeros(1,samples);
mag_drift = zeros(1,samples); ang_drift = zeros(1,samples);
noise_percentage = 0.05;

%Setup for CPU-parallellized 'parsim' function
model = 'noisy_id_simulink';
in(1:samples) = Simulink.SimulationInput(model);
load_system(model)
for i=(1:samples)
    w=W(i); dt_max = 2*pi/(w*40);
    Tstop=max(1000, 4*2*pi/w);    
    in(i) = in(i).setModelParameter('StopTime', num2str(Tstop), 'MaxStep', num2str(dt_max));
    in(i) = in(i).setVariable('w', w);
    in(i) = in(i).setVariable('noise', noise_percentage);
    in(i) = in(i).setVariable('noiseT', dt_max);
end
out = parsim(in, 'ShowProgress', 'on');


%Compute magnitude and phase angle for each frequency for true TF, TF with
%gaussian noise, and TF with drift
for i=(1:samples)
    w = W(i);
    y=out(i).y.data(:); noise=out(i).noise.data; x=out(i).x.data; t=out(i).tout;
    i_ss = floor(0.7*numel(t)); i_max = numel(t);
    y_ss = y(i_ss : i_max);
    t_ss = t(i_ss : i_max);
    
    % sin(omega.t+phi) = cos(phi)sin(omega.t) + sin(phi)cos(omega.t)
    A = [sin(w*t_ss) cos(w*t_ss) ones(size(t_ss))];
    c = A \ y_ss;
    mag(i) = norm(c(1:2));
    ang(i) = atan2(c(2), c(1));

    noise_i = randi(numel(noise),size(y));
    y_noisy = y + mag(i)*noise(noise_i);
    y_noisy_ss = y_noisy(i_ss : i_max);
    
    % sin(w.t+phi) + k = cos(phi)sin(w.t) + sin(phi)cos(w.t) + k
    A = [sin(w*t_ss) cos(w*t_ss) ones(size(t_ss))];
    c = A \ y_noisy_ss;
    mag_noisy(i) = norm(c(1:2));
    ang_noisy(i) = atan2(c(2), c(1));

    drift = t / t(end) * mag(i);
    y_drift = y + drift;
    y_drift_ss = y_drift(i_ss : i_max);

    % sin(w.t+phi) + k = cos(phi)sin(w.t) + sin(phi)cos(w.t) + k
    A = [sin(w*t_ss) cos(w*t_ss) ones(size(t_ss))];
    c = A \ y_drift_ss;
    mag_drift(i) = norm(c(1:2));
    ang_drift(i) = atan2(c(2), c(1));

    fprintf("Done: %d/%d\n",i,samples)
end
ang = unwrap(ang);
ang_noisy = unwrap(ang_noisy);
ang_drift = unwrap(ang_drift);

%Plot time response
figure
plot(t, y); hold on
plot(t, y_noisy); hold on
plot(t, y_drift); hold off
title('Time Response')
legend('Output','Noisy Output', 'Drift Output')

%Plot magnitude
figure
semilogx(W,20*log10(mag)); hold on
semilogx(W,20*log10(mag_noisy), '--'); hold on
semilogx(W,20*log10(mag_drift), '--'); hold off
title('Magnitude')
legend('Experimental', 'Experimental Noisy', 'Experimental Drift')

%Plot phase
figure
semilogx(W, ang*180/pi); hold on
semilogx(W, ang_noisy*180/pi, '--'); hold on
semilogx(W, ang_drift*180/pi, '--'); hold off
title('Phase')
legend('Experimental', 'Experimental Noisy', 'Experimental Drift')