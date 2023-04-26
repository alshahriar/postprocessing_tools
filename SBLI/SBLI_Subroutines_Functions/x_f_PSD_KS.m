function [f,PSD] = x_f_PSD_KS(signals,t)
L = length(t);
nS = length(signals(1,:));
signal = zeros(L,1);
dt = t(2)-t(1);
Fs = 1/dt;
n = 2^nextpow2(L);
f = Fs*(0:(n/2))/n;
PSD = zeros(n/2+1,nS);
for j = 1:nS
    signal(:)= signals(1:L,j);
    Y = fft(signal,n);
    P = abs(Y/n).^2;
    PSD(:,j) = P(1:n/2+1);
end
end