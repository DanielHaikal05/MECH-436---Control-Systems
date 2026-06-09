clear all; close all;

%Reference TF
W = logspace(-5, 0, 40);
num = [1 1.3893e-3 5.3055e-8];
den = [1 3.87704e-2 2.551212584e-4 8.7843039904e-8 5.619890472e-12];
K = 1.75e-4;
G=K*tf(num, den, 'InputDelay', 9.3);
[mag, ang, w_out] = bode(G, W);
mag = squeeze(mag); ang = squeeze(ang);

%Pade order 1 magnitude and phase
[num_pade1, den_pade1] = pade(9.3, 1);
G_pade1 = K * tf(num_pade1, den_pade1) * tf(num, den);
[mag_pade1, ang_pade1, w_pade1] = bode(G_pade1, W);
mag_pade1 = squeeze(mag_pade1); ang_pade1 = squeeze(ang_pade1);

%Pade order 2 magnitude and phase
[num_pade2, den_pade2] = pade(9.3, 2);
G_pade2 = K * tf(num_pade2, den_pade2) * tf(num, den);
[mag_pade2, ang_pade2, w_pade2] = bode(G_pade2, W);
mag_pade2 = squeeze(mag_pade2); ang_pade2 = squeeze(ang_pade2);

%Pade order 3 magnitude and phase
[num_pade3, den_pade3] = pade(9.3, 3);
G_pade3 = K * tf(num_pade3, den_pade3) * tf(num, den);
[mag_pade3, ang_pade3, w_pade3] = bode(G_pade3, W);
mag_pade3 = squeeze(mag_pade3); ang_pade3 = squeeze(ang_pade3);

%Pade order 4 magnitude and phase
[num_pade4, den_pade4] = pade(9.3, 4);
G_pade4 = K * tf(num_pade4, den_pade4) * tf(num, den);
[mag_pade4, ang_pade4, w_pade4] = bode(G_pade4, W);
mag_pade4 = squeeze(mag_pade4); ang_pade4 = squeeze(ang_pade4);


%Magnitude plot
figure
semilogx(w_out,20*log10(mag)); hold on
semilogx(w_pade1,20*log10(mag_pade1), '--'); hold on
semilogx(w_pade2, 20*log10(mag_pade2), '--'); hold on
semilogx(w_pade3, 20*log10(mag_pade3), '--'); hold on
semilogx(w_pade4, 20*log10(mag_pade4), '--'); hold off
title('Magnitude')
legend('G', 'Pade Order 1', 'Pade Order 2', 'Pade Order 3', 'Pade Order 4')

%Phase plot
figure
semilogx(w_out,ang); hold on
semilogx(w_pade1,ang_pade1-360, '--'); hold on
semilogx(w_pade2, ang_pade2-360, '--'); hold on
semilogx(w_pade3, ang_pade3-720, '--'); hold on
semilogx(w_pade4, ang_pade4-720, '--'); hold off
title('Phase')
legend('G', 'Pade Order 1', 'Pade Order 2', 'Pade Order 3', 'Pade Order 4')

%Root Loci
figure; rlocus(G_pade1); title('Order 1')
figure; rlocus(G_pade2); title('Order 2')
figure; rlocus(G_pade3); title('Order 3')
figure; rlocus(G_pade4); title('Order 4')


%Linearize simulink model
[A,B,C,D] = linmod('linearize_simulink');
[num_lin, den_lin] = ss2tf(A,B,C,D);
G_lin = tf(num_lin, den_lin);

p = eig(A);
z = roots(num_lin);