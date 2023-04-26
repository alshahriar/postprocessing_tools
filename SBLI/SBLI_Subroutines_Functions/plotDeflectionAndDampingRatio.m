% seperate the signal in different parts
% calculating the frequency and damping ration for each part
function [peaks,peaksTimeIndex,valleys,valleysTimeIndex,dampingRatio,freqValues,figNo] = plotDeflectionAndDampingRatio(caseName,imageFormat,defAtMiddle,t_norm,nParts)
eachPart = floor(length(defAtMiddle)/nParts);
for i=1:nParts
    startPointx(i) = 1+eachPart*(i-1);
    endPointx(i) = eachPart*i;
end
endPointx(end) = length(defAtMiddle);
countp = 0;
countv = 0;
countIndexP = 0;
countIndexV = 0;

defAtMiddleNorm = defAtMiddle - mean(defAtMiddle);
%defAtMiddleNorm = smooth(defAtMiddleNorm,15);
% plotting def at mid
figNo=1234564;
figure(figNo); clf;
for i = 1:nParts
    t_norm2 = t_norm(startPointx(i):endPointx(i));
    defAtMiddle2 = (defAtMiddleNorm(startPointx(i):endPointx(i)));
    tempAverageDefatMiddle = mean(defAtMiddle2);
    maxDef = max(defAtMiddle2);
    %minDef = min(defAtMiddle);
    [peaksTemp,peaksTimeIndexTemp] = findpeaks(defAtMiddle2,'MinPeakHeight',tempAverageDefatMiddle);
    invertedDef = maxDef - defAtMiddle2;
    [valleysTemp,valleysTimeIndexTemp] = findpeaks(invertedDef,'MinPeakHeight',maxDef - tempAverageDefatMiddle);
    valleysTemp = maxDef - valleysTemp;
    nPeaksTemp = length(peaksTemp);
    nValleysTemp = length(valleysTemp);
    peaks(countp+1:countp+nPeaksTemp) = peaksTemp;
    valleys(countv+1:countv+nValleysTemp) = valleysTemp;
    peaksTimeIndex(countp+1:countp+nPeaksTemp) = countIndexP + peaksTimeIndexTemp;
    valleysTimeIndex(countv+1:countv+nValleysTemp) = countIndexV + valleysTimeIndexTemp;
    countIndexP = countIndexP + length(defAtMiddle2);
    countIndexV = countIndexV + length(defAtMiddle2);
    countp = countp + nPeaksTemp;
    countv = countv + nValleysTemp;
    subplot(3,1,1)
    %     plot(t_norm2,defAtMiddle2)
    hold on;
    plot(t_norm2,defAtMiddleNorm(startPointx(i):endPointx(i)),'-k')
    plot(t_norm2(peaksTimeIndexTemp),peaksTemp,'*r');
    plot(t_norm2(valleysTimeIndexTemp),valleysTemp,'*b')

    clear t_norm2 defAtMiddle2
end
    subplot(3,1,1)
    hold on



    title("deflection")

clear indexLocation dampingRatio
nPeaks = length(peaksTimeIndex);
x1 = defAtMiddleNorm(peaksTimeIndex(end));
for i = 1:nPeaks-1
    x0 = defAtMiddleNorm(peaksTimeIndex(i));
    nPeaksForDamping = nPeaks - i;
    dampingRatio(1,i) = t_norm(peaksTimeIndex(i));
    dampingRatio(2,i) = 1/nPeaksForDamping * log(x0/x1);
end

nPeaksFFT = 3; % number of peaks for FFT
nPeaksFFT = min(nPeaksFFT,length(peaks)-1);
for i = 1: length(peaks)-nPeaksFFT
    signalFFT2 = defAtMiddleNorm(peaksTimeIndex(i):peaksTimeIndex(i+nPeaksFFT));
    t_FFT = t_norm(peaksTimeIndex(i):peaksTimeIndex(i+nPeaksFFT));
    [f_2,xt_2] = x_f_PSD_2(t_FFT,signalFFT2);
    [~,IndexTemp] = max(xt_2);
    freqValues(2,i) = f_2(IndexTemp);
    freqValues(1,i) = t_norm(peaksTimeIndex(i));
end

for i = 1: length(peaks)
    x0 = defAtMiddleNorm(peaksTimeIndex(i));
    hold on
    %     text(t_norm(peaksTimeIndex(i)),x0,'>')
end

for i = 1: length(valleys)
    x0 = defAtMiddleNorm(valleysTimeIndex(i));
    hold on
    %     text(t_norm(valleysTimeIndex(i)),x0,'<')
end

% Plotting damping ratio
subplot(3,1,2)
plot(dampingRatio(1,:),dampingRatio(2,:))
title("damping ratio")
% plotting frequency - kind of specturm graph
subplot(3,1,3)
plot(freqValues(1,:),freqValues(2,:))
title("frequency")
end
