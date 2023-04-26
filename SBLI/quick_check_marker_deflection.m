clear all

if(isunix)
    addpath('SBLI_Subroutines_Functions')
end

imageFormat = 'png';
% imageFormat = 'epsc';

if(0)
    % normalPanel
    xMarkerS = -14.8;
    xMarkerE = 14.1;
    xMarkerFlexS = -5;
    xMarkerFlexE = 5;
elseif(0)
    %longerPanel
    % xMarkerS = -19;
    % xMarkerE = 19;
    % xMarkerFlexS = -7.5;
    % xMarkerFlexE = 7.5;
elseif(1)
    % Visbal setup
    xMarkerS = -0.5;
    xMarkerE = 2.1;
    xMarkerFlexS = 0;
    xMarkerFlexE = 1;
end

maxDeflectionValue = 0.5;

gamma = 1.4;
pinf = 1/gamma;
rhoinf = 1;
timeIndexForRigid = 75+127+183;
T_r = 1.678;
panelThickness = 0.5;
panelLenght = 10;
panelDensity = 1000;        
panelLocationY = 1;
poissionRatio = 0.45;

n_T = [0.6 1 1.34 1.678]; % for adiabatic - n_T = 0;
n_E = [5000 10000 50000];
if(isunix)
    dirNameParametricStudy = '/gpfs/research/engineering/as17r/ParametricStudy/AllCases2/latestVersion';
else
    dirNameParametricStudy = 'D:\CCNS\parametricStudy\AllCases2';
end
subfolder = "flow_output";

w0_p = (pi^2) * sqrt(n_E(1)*panelThickness^2/(12*(1-poissionRatio^2)*panelDensity*panelLenght^4)); % natural frequency
w0_c = 22.4*sqrt((n_E(1) * panelThickness^2)/(12*(1-poissionRatio^2)*panelDensity*panelLenght^4))/(2*pi); %Roark's Formulas for Stress and Strain

only_user_defined_dir = 1;
dirAppendFlag = 0;
userDirCaseNumber = 0;
if(only_user_defined_dir==1)
    userDirCaseNumber = userDirCaseNumber + 1;
    if(isunix)
        uDirName{userDirCaseNumber} = '/gpfs/research/engineering/as17r/ParametricStudy/AllCases2/latestVersion';
    else
        uDirName{userDirCaseNumber} = 'D:\CCNS\parametricStudy\AllCases2';
    end
    uNameOfFolders{userDirCaseNumber} = 'newTahoeRigid28f';
    un_T(userDirCaseNumber) = 0.6;
    un_E(userDirCaseNumber) = 5000;
else
    uDirName = [];
    uNameOfFolders = [];
    un_T = [];
    un_E = [];
end
unT = length(un_T);
unE = length(un_T);

[dirName,folderNames,n_T_Case,~] = getFileAndFolderName(n_T,n_E,dirNameParametricStudy,only_user_defined_dir,dirAppendFlag,uDirName,uNameOfFolders,un_T,un_E);

nFolder = length(folderNames);
totalNumberOfCases = length(folderNames);
%%
% prompt = ("what is the case number[1-15]?");
% caseNumber = input(prompt);
for caseNumber = 1:totalNumberOfCases
    clear defAtMiddle xLoc pAtMiddle
    load(fullfile(dirName{caseNumber},folderNames{caseNumber},subfolder,'marker.mat'),'marker');
    
    [simParameters] = loadCaseData(dirName{caseNumber},folderNames{caseNumber});
    [tahoeParam] = loadTahoeInput(dirName{caseNumber},folderNames{caseNumber});
    
    delta = simParameters.delta;
    CauchyNumber(caseNumber) = rhoinf*(simParameters.Uinf^2)/tahoeParam.E;
    if(1)
        originalMarker = marker; clear marker;
        % Removing the part that is not flexibile
        marker = originalMarker(:,:,timeIndexForRigid+1:end);
    end
    nTime = length(marker(1,1,:));
    t_norm = simParameters.ntec*simParameters.dt*simParameters.Uinf/delta*(1:nTime);
    t_norm = simParameters.ntec*simParameters.dt*(1:nTime);

    xx = marker(:,1,1);
    midLocation_x = 0.5*(xMarkerFlexE+xMarkerFlexS);
    three4th_Location_x = xMarkerFlexS + 0.75*(xMarkerFlexE-xMarkerFlexS);
    [~,midPointIndex] = min(abs(xx-midLocation_x));
    [~,three4thIndex] = min(abs(xx-three4th_Location_x));
    disp(['Case: ',num2str(caseNumber),', Name: ', folderNames{caseNumber},', x-coord of midpoint is:', num2str(marker(midPointIndex,1,1))])
    defAtMiddle(:) = marker(midPointIndex,2,:);
    defAt3_4th(:) = marker(three4thIndex,2,:);    
    pAtMiddle(:) = marker(midPointIndex,8,:);
    
    xLoc(:) = marker(midPointIndex,1,:);
    averageDefatMiddle = mean(defAtMiddle(1:end));
    maxDef = max(max(marker(:,2,1:end)));
    minDef = min(min(marker(:,2,1:end)));
    threshold = (maxDef - minDef) * 0.03;
    [peaks,peaksTimeIndex] = findpeaks(defAtMiddle,'MinPeakHeight',averageDefatMiddle);
    invertedDef = maxDef - defAtMiddle;
    [valleys,valleysTimeIndex] = findpeaks(invertedDef,'MinPeakHeight',maxDef - averageDefatMiddle);
    valleys = maxDef - valleys;
    
    J = customcolormap([0 0.5 1], {'#00FF00','#ff0000','#0000FF'},nTime);   %
    color = colormap(J);
    color = [color];
    
    rowsPlot = 2;
    colPlot = 2;
    figure; clf;
    subplot(rowsPlot,colPlot,1);
    step = 2;
    for i = 1:step:nTime-step
        plot(t_norm(i:i+step),defAtMiddle(i:i+step),'Color',color(i,:));
        hold on
    end
    for i = 1:step:nTime-step
        plot(t_norm(i:i+step),defAt3_4th(i:i+step),':','Color',color(i,:));
        hold on
    end
    % text(dt*peaksTimeIndex,peaks,'*')
    % text(dt*valleysTimeIndex,valleys,'*')
    titleText = sprintf('T_w = %2.2fT_r and E = %d (Case:%d)',n_T_Case(caseNumber),tahoeParam.E,caseNumber);
    title(['Oscillation at the Middle of the Panel ',titleText])
    xlabel('time, t/')
    ylabel('y/')
    set(gca, 'FontName', 'Times')
    
    subplot(rowsPlot,colPlot,2);
    step = 2;
    for i = 1:step:nTime-step
        plot(xLoc(i:i+step),defAtMiddle(i:i+step),'Color',color(i,:));
        hold on
    end
    % text(dt*peaksTimeIndex,peaks,'*')
    % text(dt*valleysTimeIndex,valleys,'*')
    titleText = sprintf('T_w = %2.2fT_r and E = %d (Case:%d)',n_T_Case(caseNumber),tahoeParam.E,caseNumber);
    title(['movement of the mid point',titleText])
    xlabel('x/')
    ylabel('y/')
    set(gca, 'FontName', 'Times')
    
    % Deflection profile at the 2nd vallays (max delfection) and 2nd peaks(min deflection)
    subplot(rowsPlot,colPlot,3);
    for i=1:length(peaks)
        plot(marker(:,1,peaksTimeIndex(i)), marker(:,2,peaksTimeIndex(i)),'Color',color(peaksTimeIndex(i),:))
        hold on
    end
    title("Panel Profile at Min Deflection")
    xlabel('x/')
    ylabel('y/')
    axis([xMarkerFlexS xMarkerFlexE -inf inf])
    set(gca, 'FontName', 'Times')
    
    subplot(rowsPlot,colPlot,4);
    for i=1:length(valleys)
        plot(marker(:,1,valleysTimeIndex(i)), marker(:,2,valleysTimeIndex(i)),'Color',color(valleysTimeIndex(i),:))
        hold on
    end
    title("Panel Profile at Max Deflection")
    xlabel('x/')
    ylabel('y/')
    axis([xMarkerFlexS xMarkerFlexE minDef maxDef])
    set(gca, 'FontName', 'Times')
    set(gcf, 'Position',  [100, 100, 800, 300])
    if(isunix)
        saveas(gcf,"n_SpecificCases_"+caseNumber+"_"+folderNames{caseNumber},'epsc')
    end
    defAtMiddleAllCases{caseNumber}(:,1) = t_norm;
    defAtMiddleAllCases{caseNumber}(:,2) = defAtMiddle;
    defAtMiddleAllCases{caseNumber}(:,3) = pAtMiddle;   
end

%% In frequency domain
if(1)
    figure;
    for caseNumber = 1:totalNumberOfCases
        [simParameters] = loadCaseData(dirName{caseNumber},folderNames{caseNumber});
        t = defAtMiddleAllCases{caseNumber}(:,1);
        d = defAtMiddleAllCases{caseNumber}(:,2);
        [freqvec,xdft] = x_f_PSD_2(t,d);
        plot(freqvec,abs(xdft)); hold on
        [sigr, coeffiecients] = fit_damped_sinewave((d-mean(d))');
        coeffiecients(2) = simParameters.ntec*simParameters.dt*coeffiecients(2)*delta/simParameters.Uinf;
        coeffiecients(3) = simParameters.ntec*simParameters.dt*coeffiecients(3)*delta/simParameters.Uinf;
        dampingRatio2 = coeffiecients(2)/sqrt(coeffiecients(2)^2+coeffiecients(3)^2);
        sprintf('A = %f, alpha = %f, w/(2pi) = %f, phi = %f, xi = %f',coeffiecients,dampingRatio2)
    end
    xlim([0 .0025]);
    ylim([0 12]);
    legend('$P_c/P_\infty$ = 1.0','$P_c/P_\infty$ = 1.2','$P_c/P_\infty$ = 1.4','$P_c/P_\infty = \langle P \rangle$','Interpreter','latex')
    set(gca,'FontName','Times')
end

%%
if(totalNumberOfCases==3)
    colorS = ['r';'b';'g'];
elseif(totalNumberOfCases==4)
    colorS = ['r';'b';'g';'k'];
elseif(totalNumberOfCases==5)
    colorS = ['r';'b';'g';'k';'m'];
elseif(totalNumberOfCases==2)
    colorS = ['r';'b'];
else
    colorS = ['k'];
end

figure;
for caseNumber = 1:totalNumberOfCases
    [simParameters] = loadCaseData(dirName{caseNumber},folderNames{caseNumber});
    t = defAtMiddleAllCases{caseNumber}(1:end,1);
    d = defAtMiddleAllCases{caseNumber}(1:end,3);
    plot(t,(d-panelLocationY)/delta,colorS(caseNumber));
    hold on;
end
% xlim([200 20960]);
% legend('\rho_s = 1000','\rho_s = 500','\rho_s = 100')
legend('$P_c/P_\infty$ = 1.0','$P_c/P_\infty$ = 1.2','$P_c/P_\infty$ = 1.4','$P_c/P_\infty = \langle P \rangle$','Interpreter','latex')
xlabel('time, $tU/\delta$','Interpreter','latex')
ylabel('$y/\delta$','Interpreter','latex')
set(gca, 'FontName', 'Times')
set(gcf, 'Position',  [100, 100, 400, 300])
grid on

%% FFT at 3 different time
if(0)
    caseNumber = 1;
    [simParameters] = loadCaseData(dirName{caseNumber},folderNames{caseNumber});
    t = defAtMiddleAllCases{caseNumber}(:,1);
    d = defAtMiddleAllCases{caseNumber}(:,2);
    
    % Part - 1
    t1v = 15480;
    t2v = 36120;
    
    sIndex = 1;
    eIndex = find(t==t1v);
    
    tt = t(sIndex:eIndex);
    dd = d(sIndex:eIndex);
    [freqvec,xdft] = x_f_PSD_2(tt,dd);
    figure;
    plot(freqvec,abs(xdft),'-r')
    title('Part 1');
    
    
    % Part - 2
    sIndex = eIndex+1;
    eIndex = find(t==t2v);
    tt = t(sIndex:eIndex);
    dd = d(sIndex:eIndex);
    [freqvec,xdft] = x_f_PSD_2(tt,dd);
    figure;
    plot(freqvec,abs(xdft),'-r')
    title('Part 2');
    
    % Part - 3
    sIndex = eIndex+1;
    eIndex = length(t);
    tt = t(sIndex:eIndex);
    dd = d(sIndex:eIndex);
    [freqvec,xdft] = x_f_PSD_2(tt,dd);
    figure;
    plot(freqvec,abs(xdft),'-r')
    title('Part 3');
end
