function [f,xdft] = x_f_PSD(tt,tsignals)
n = length(tt);
if(mod(n,2)~=0)
    t = tt(1:n-1);
    signals(:,:) = tsignals(:,1:n-1);
    n = n - 1;
else
    t = tt;
    signals = tsignals;
end

signal = zeros(n,1);
nS = length(signals(:,1));

dt = t(2)-t(1);
fs = 1/dt;
fshift = (-n/2:n/2-1)*(fs/n);
[V,I] = max(fshift==0);

f = fshift(I:n);
xdft = zeros(nS,n-I+1);

for j = 1:nS
    for i = 1:n
        signal(i)= signals(j,i);
    end
    ytemp = fft(signal);
    yshift = fftshift(ytemp);
    y = abs(yshift(I:n));
    y(1) = 0;
    xdft(j,:) = y;
end
end