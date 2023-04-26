function [f,y] = x_f_PSD_2(tt,tsignal)
n = length(tt);
if(mod(n,2)~=0)
    t = tt(1:n-1);
    signal(:) = tsignal(1:n-1);
    n = n - 1;
else
    t = tt;
    signal = tsignal;
end

dt = t(2)-t(1);
fs = 1/dt;
fshift = (-n/2:n/2-1)*(fs/n);
[V,I] = max(fshift==0);
f = fshift(I:n);

ytemp = fft(signal);   
yshift = fftshift(ytemp);
y = abs(yshift(I:n));
y(1) = y(2);
end



