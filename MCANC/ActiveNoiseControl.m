% A simplified version of single-channel and multiple-channel ANC systems 

close all; clear; clc;

T = 1e3;  % number of points
L = 16;  % order of the filter
Fs = 5e3;  % sampling frequency
t = (0:T-1) / Fs;  % time in second
S = 5e6 * sin(2 * pi * 113 * t) + 4.78e6 * sin(2 * pi * 228 * t) + 5.12e6 * sin(2 * pi * 311 * t);  % sine wave
L1 = 2;  % length in meter
L2 = 2;
r1 = 4;
r2 = 4;
phi1 = 90;  % angle in degree
phi2 = 0;
the1 = 60;
the2 = 20;
c = 345;  % sound speed(m/s)
r11 = sqrt((r1 * cosd(the1) - L1 * cosd(phi1))^2 + (r1 * sind(the1) - L1 * sind(phi1))^2);  % y1 to e1
r22 = sqrt((r2 * cosd(the2) - L2 * cosd(phi2))^2 + (r2 * sind(the2) - L2 * sind(phi2))^2);  % y2 to e2
r21 = sqrt((r1 * cosd(the1) - L2 * cosd(phi2))^2 + (r1 * sind(the1) - L2 * sind(phi2))^2);  % y2 to e1
r12 = sqrt((r2 * cosd(the2) - L1 * cosd(phi1))^2 + (r2 * sind(the2) - L1 * sind(phi1))^2);  % y1 to e2
tr1 = r1 / c;  % required time for traveling
tr2 = r2 / c;
tr11 = r11 / c;
tr22 = r22 / c;
tr12 = r12 / c;
tr21 = r21 / c;
SA = sqrt((r1 * cosd(the1) + 0.5 * (r2 * cosd(the2) - r1 * cosd(the1)))^2 + (r1 * sind(the1))^2);  % distance in meter
SB = sqrt((r2 * cosd(the2))^2 + (r1 * sind(the1))^2);
SC = sqrt((r2 * cosd(the2))^2 + (r2 * sind(the2) + 0.5 * (r1 * sind(the1) - r2 * sind(the2)))^2);
SD = sqrt((r1 * cosd(the1) + 0.5 * (r2 * cosd(the2) - r1 * cosd(the1)))^2 + (r2 * sind(the2))^2);
SE = sqrt((r1 * cosd(the1))^2 + (r2 * sind(the2))^2);
SF = sqrt((r1 * cosd(the1))^2 + (r2 * sind(the2) + 0.5 * (r1 * sind(the1) - r2 * sind(the2)))^2);
Y1A = sqrt((r1 * cosd(the1) + 0.5 * (r2 * cosd(the2) - r1 * cosd(the1)) - L1 * cosd(phi1))^2 + (r1 * sind(the1) - L1 * sind(phi1))^2);
Y1B = sqrt((r2 * cosd(the2) - L1 * cosd(phi1))^2 + (r1 * sind(the1) - L1 * sind(phi1))^2);
Y1C = sqrt((r2 * cosd(the2) - L1 * cosd(phi1))^2 + (r2 * sind(the2) + 0.5 * (r1 * sind(the1) - r2 * sind(the2)) - L1 * sind(phi1))^2);
Y1D = sqrt((r1 * cosd(the1) + 0.5 * (r2 * cosd(the2) - r1 * cosd(the1)) - L1 * cosd(phi1))^2 + (r2 * sind(the2) - L1 * sind(phi1))^2);
Y1E = sqrt((r1 * cosd(the1) - L1 * cosd(phi1))^2 + (r2 * sind(the2) - L1 * sind(phi1))^2);
Y1F = sqrt((r1 * cosd(the1) - L1 * cosd(phi1))^2 + (r2 * sind(the2) + 0.5 * (r1 * sind(the1) - r2 * sind(the2)) - L1 * sind(phi1))^2);
Y2A = sqrt((r1 * cosd(the1) + 0.5 * (r2 * cosd(the2) - r1 * cosd(the1)) - L2 * cosd(phi2))^2 + (r1 * sind(the1) - L2 * sind(phi2))^2);
Y2B = sqrt((r2 * cosd(the2) - L2 * cosd(phi2))^2 + (r1 * sind(the1) - L2 * sind(phi2))^2);
Y2C = sqrt((r2 * cosd(the2) - L2 * cosd(phi2))^2 + (r2 * sind(the2) + 0.5 * (r1 * sind(the1) - r2 * sind(the2)) - L2 * sind(phi2))^2);
Y2D = sqrt((r1 * cosd(the1) + 0.5 * (r2 * cosd(the2) - r1 * cosd(the1)) - L2 * cosd(phi2))^2 + (r2 * sind(the2) - L2 * sind(phi2))^2);
Y2E = sqrt((r1 * cosd(the1) - L2 * cosd(phi2))^2 + (r2 * sind(the2) - L2 * sind(phi2))^2);
Y2F = sqrt((r1 * cosd(the1) - L2 * cosd(phi2))^2 + (r2 * sind(the2) + 0.5 * (r1 * sind(the1) - r2 * sind(the2)) - L2 * sind(phi2))^2);
E1E = sqrt((r2 * sind(the2) - r1 * sind(the1))^2);
E1B = sqrt((r2 * cosd(the2) - r1 * cosd(the1))^2);
A = E1E * E1B;

% secondary path offline training
% white noise signal
xs=wgn(1, T, 0);

% measure the reference signal, delay is because of the distance between them
ds11 = delayseq(xs / r11, round(tr11), Fs);
ds21 = delayseq(xs / r21, round(tr21), Fs);
ds12 = delayseq(xs / r12, round(tr12), Fs);
ds22 = delayseq(xs / r22, round(tr22), Fs);

% initiate the system
Shx = zeros(1, L);
% the weight of Sh 
Sh11 = Shx;          
Sh21 = Shx;  % y1 to e2
Sh12 = Shx;  % y2 to e1
Sh22 = Shx;
es11 = zeros(1, T);  % error
es21 = es11;
es12 = es11;
es22 = es11;

mu1 = 0.01;  % learning rate
for n = 1:T
    Shx = [xs(n) Shx(1:L-1)];  % update the state
    es11(n) = ds11(n) - sum(Shx .* Sh11);  % calculate the error
    es21(n) = ds21(n) - sum(Shx .* Sh21);
    es12(n) = ds12(n) - sum(Shx .* Sh12);
    es22(n) = ds22(n) - sum(Shx .* Sh22);
    Sh11 = Sh11 + mu1 * es11(n) * Shx;  % adjust the weight of Sh
    Sh21 = Sh21 + mu1 * es21(n) * Shx;
    Sh12 = Sh12 + mu1 * es12(n) * Shx;
    Sh22 = Sh22 + mu1 * es22(n) * Shx;
end

% ANC
snr = 10;  % SNR
X = awgn(S, snr, 'measured');  % signal w background noise

% measure the reference signal
d1 = delayseq(S / r1, round(tr1), Fs);
d2 = delayseq(S / r2, round(tr2), Fs);
d = d1;

% initiate the system
x = zeros(1, L);  % input

% weight
% 1x1x1
w = zeros(1, L);
% 1x2x2
w1 = w;
w2 = w;

% speaker
% 1x1x1
y = zeros(1, T);
% 1x2x2
y1 = y;
y2 = y;

% error mic
% 1x1x1
e = zeros(1, T);
% 1x2x2
e1 = e;
e2 = e;

% filtered x
% 1x1x1
xp = zeros(1, L);
% 1x2x2
xp11 = xp;
xp12 = xp;
xp21 = xp;
xp22 = xp;

mu2 = 1e-15;  % learning rate
% ANC 1x1x1
for n = 1:T  % discrete time n
    x = [X(n) x(1:L-1)];  % update the controller state    
    y(n) = sum(x .* w);  % propagate to secondary path
    e(n) = awgn(d(n) - delayseq(y(n) / r11, round(tr11), Fs), snr, 'measured');  % measure the residue
    xp = [sum(x .* Sh11) xp(1:L-1)];  % update filtered x
    w = w + mu2 * xp * e(n);  % adjust the weight
end

% ANC 1x2x2
for n = 1:T  % discrete time n
    x = [X(n) x(1:L-1)];  % update the controller state    
    y1(n) = sum(x .* w1);  % propagate to secondary path
    y2(n) = sum(x .* w2);
    e1(n) = awgn(d1(n) - delayseq(y1(n) / r11, round(tr11), Fs) - delayseq(y2(n) / r21, round(tr21), Fs), snr, 'measured');  % measure the residue
    e2(n) = awgn(d2(n) - delayseq(y2(n) / r22, round(tr22), Fs) - delayseq(y1(n) / r12, round(tr12), Fs), snr, 'measured');
    xp11 = [sum(x .* Sh11) xp11(1:L-1)];  % update filtered x
    xp12 = [sum(x .* Sh21) xp12(1:L-1)];
    xp21 = [sum(x .* Sh12) xp21(1:L-1)];
    xp22 = [sum(x .* Sh22) xp22(1:L-1)];
    w1 = w1 + mu2 * e1(n) * (1 * xp11 + 0 * xp12);  % adjust the weight
    w2 = w2 + mu2 * e2(n) * (1 * xp22 + 0 * xp21);
end

% plot
figure('name', 'Multiple ANC 1x2x2')
plot([1:T], pow2db(abs(X)), 'k', [1:T], pow2db(abs(e1)))
title('1x2x2')
ylabel('dB');
xlabel('Discrete time n');
legend('Noise input', 'Noise residue')
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 20)

figure('name', 'ANC 1x1x1')
plot([1:T], pow2db(abs(X)), 'k', [1:T], pow2db(abs(e)))
title('1x1x1')
ylabel('dB');
xlabel('Discrete time n');
legend('Noise input', 'Noise residue')
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 20)

figure('name', 'Juxtaposition')
subplot(3, 1, 1)
plot([1:T], e)
ylabel('Amplitude');
xlabel('Discrete time n');
legend('Noise residue(1x1x1)')
subplot(3, 1, 2)
plot([1:T], e1)
ylabel('Amplitude');
xlabel('Discrete time n');
legend('Noise residue(1x2x2)')
subplot(3, 1, 3 )
plot([1:T], X)
ylabel('Amplitude');
xlabel('Discrete time n');
legend('Noise input')
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 15)

figure('name','FFT')
Xfft = fft(X);
e1fft = fft(e1);
XP2 = abs(Xfft/T);
XP1 = XP2(1:T/2+1);
XP1(2:end-1) = 2 * XP1(2:end-1);
e1P2 = abs(e1fft / T);
e1P1 = e1P2(1:T / 2+1);
e1P1(2:end-1) = 2 * e1P1(2:end-1);
f = Fs * (0:(T/2)) / T;
plot(f, pow2db(XP1), f, pow2db(e1P1))
title('FFT(1x2x2)')
xlabel('f(Hz)')
ylabel('dB')
legend('Noise input','Noise residue')
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 20)

figure('name', 'FFT')
efft = fft(e);
eP2 = abs(efft/T);
eP1 = eP2(1:T/2+1);
eP1(2:end-1) = 2 * eP1(2:end-1);
plot(f, pow2db(XP1), f, pow2db(eP1))
title('FFT(1x1x1)')
xlabel('f(Hz)')
ylabel('dB')
legend('Noise input', 'Noise residue')
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 20)