function [f, P1] = fft_P1(X, Fs)
% Nyquist Amplitude of X

% Fs, Nyquist sampling Frequence
% X, signal.

Y = fft(X);
L = length(X);
P2 = abs(Y/L);
P1 = P2(1:ceil(L/2+1));
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:ceil(L/2))/L;
% plot(f, P1)
% xlim([0.1 f(end)])
% title('Single-Sided Amplitude Spectrum of X(t)')
% xlabel('f (Hz)')
% ylabel('|P1(f)|')