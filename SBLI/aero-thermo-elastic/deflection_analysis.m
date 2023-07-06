%% Section 1


clear all

goAhead = input('Are you sure?');

addpath('/gpfs/research/engineering/as17r/postprocessing_tools/SBLI/SBLI_Subroutines_Functions');

imageFormat = 'png';
% imageFormat = 'epsc';

if(1)
    % normalPanel
    xMarkerS = -14.8;
    xMarkerE = 14.1;
    xMarkerFlexS = -4.98;
    xMarkerFlexE = 4.98;
    panelLocationY = 1;
    timeIndexForRigid = 0;

    panelThickness = 0.5;
    panelLenght = xMarkerFlexE-xMarkerFlexS;

elseif(0)
    %longerPanel
    % xMarkerS = -19;
    % xMarkerE = 19;
    % xMarkerFlexS = -7.5;
    % xMarkerFlexE = 7.5;
elseif(0)
    % Visbal setup
    xMarkerS = -0.5;
    xMarkerE = 2.1;
    xMarkerFlexS = 0;
    xMarkerFlexE = 1;
    panelLocationY = 0;
    timeIndexForRigid = 00;
    sDefCurve = 300;
    eDefCurve = 709;
end

%longerPanel
maxDeflectionValue = 0.5;

gamma = 1.4;
pinf = 1/gamma;
rhoinf = 1;

T_r = 1.678;
userDirCaseNumber = 0;
if(isunix)

    if(0)
        case_list_old_study_set_1;
        case_list_thin_panel_set_1;
        case_list_thermo_set_1;
        case_list_thermo_set_2;
        case_list_thermo_set_3;
    else
        userDirCaseNumber = userDirCaseNumber + 1;
        dirName{userDirCaseNumber} = '/gpfs/home/as17r/engineering/ParametricStudy/AllCases2/Thermal_Coupling_Test_Cases/OldSetup/';
        folderNames{userDirCaseNumber} = 'test34_2_rerun';
        iterS(userDirCaseNumber) = 40; iterE(userDirCaseNumber) = 0;
        
        userDirCaseNumber = userDirCaseNumber + 1;
        dirName{userDirCaseNumber} = '/gpfs/home/as17r/engineering/ParametricStudy/AllCases2/Thermal_Coupling_Test_Cases/OldSetup/';
        folderNames{userDirCaseNumber} = 'test33_3_rerun';
        iterS(userDirCaseNumber) = 40; iterE(userDirCaseNumber) = 0;

        for i =1:userDirCaseNumber
            uSubfolders{i} = 'flow_output';
            nvarList(i) = 11;
        end
        srow = 3056; erow = 4048;
    end



else
    dirName{userDirCaseNumber} = 'D:\CCNS\parametricStudy\AllCases2';
    folderNames{userDirCaseNumber} = 'H_T0600_E05000_veryLongRun';
    n_T_Case(userDirCaseNumber) = 0.6;
    n_E_Case(userDirCaseNumber) = 5000;

    userDirCaseNumber = userDirCaseNumber + 1;
    dirName{userDirCaseNumber} = 'D:\CCNS\parametricStudy\AllCases2';
    folderNames{userDirCaseNumber} = 'H_T0600_E05000_veryLongRun_density_100';
    n_T_Case(userDirCaseNumber) = 0.6;
    n_E_Case(userDirCaseNumber) = 5000;

    userDirCaseNumber = userDirCaseNumber + 1;
    dirName{userDirCaseNumber} = 'D:\CCNS\parametricStudy\AllCases2';
    folderNames{userDirCaseNumber} = 'H_T0600_E05000_veryLongRun_density_500_oldcode';
    n_T_Case(userDirCaseNumber) = 0.6;
    n_E_Case(userDirCaseNumber) = 5000;

end

nFolder = length(folderNames);
totalNumberOfCases = length(folderNames);

defAtMiddleAllCases = cell(totalNumberOfCases,1);
def14thAllCases = cell(totalNumberOfCases,1);
def34thAllCases = cell(totalNumberOfCases,1);
freqValuesAllCases = cell(totalNumberOfCases,1);
peaksTimeIndexAllCases = cell(totalNumberOfCases,1);
valleysTimeIndexAllCases = cell(totalNumberOfCases,1);
coeffiecientsAllCases = zeros(totalNumberOfCases,5);
t_normAllCases = cell(totalNumberOfCases,1);
dampingRatioAllCases = cell(totalNumberOfCases,1);
load(fullfile(dirName{1},folderNames{1},'marker.mat'));
multipleCasesPr = zeros(length(marker(:,1,1)),11,totalNumberOfCases);
clear tempData_1;
fname = 'dampedSinosoid_user.csv';
fidDampedSinosoid = fopen(fname,'W');
fprintf(fidDampedSinosoid,'%s,%s,%s,%s,%s,%s\n','Name','A','lambda','omega','phi','zeta');

%%
fname = 'calculatedParam.csv';
fidLambdaPaper = fopen(fname,'W');
for caseNumber = 1:totalNumberOfCases
    [simParameters] = loadCaseData(dirName{caseNumber},folderNames{caseNumber});
    [solid] = loadTahoeInput(dirName{caseNumber},folderNames{caseNumber});

    D_paper = solid.E*(panelThickness^3)/(12*(1-solid.nu^2));
    lambda_paper = rhoinf*(simParameters.Uinf^2)*(panelLenght^3)/D_paper;
    fprintf(fidLambdaPaper,'lambda = %f\n',lambda_paper);
end
fclose(fidLambdaPaper);

%%
averageDefatMiddle = zeros(1,totalNumberOfCases);
CauchyNumber = zeros(1,totalNumberOfCases);
freqValue = zeros(1,totalNumberOfCases);
%%
for caseNumber = 1:totalNumberOfCases
    caseName = folderNames{caseNumber};
    load(fullfile(dirName{caseNumber},folderNames{caseNumber},'marker'));

    xx = marker(:,1,1);
    midLocation_x = 0.5*(xMarkerFlexE+xMarkerFlexS);
    three4th_Location_x = xMarkerFlexS + 0.75*(xMarkerFlexE-xMarkerFlexS);
    one4th_Location_x = xMarkerFlexS + 0.25*(xMarkerFlexE-xMarkerFlexS);

    [~,midPointIndex] = min(abs(xx-midLocation_x));
    [~,three4thIndex] = min(abs(xx-three4th_Location_x));
    [~,one4thIndex] = min(abs(xx-one4th_Location_x));

    [fluid] = loadCaseData(dirName{caseNumber},folderNames{caseNumber});
    [solid] = loadTahoeInput(dirName{caseNumber},folderNames{caseNumber});

    if(1)
        delta = 10;
    else
        delta = fluid.delta;
    end
    disp("Check delta: "+delta);

    CauchyNumber(caseNumber) = rhoinf*(fluid.Uinf^2)/solid.E;

    %     disp(['loading: ',folderNames{caseNumber}])
    disp(['case: ',num2str(caseNumber),', Name: ', folderNames{caseNumber},', x-coord of midpoint is:', num2str(marker(midPointIndex,1,1))])
    if(1)
        figure;
        % quick check on the time series of deflection
        plot(1:length(marker(midPointIndex,2,:)),squeeze(marker(midPointIndex,2,:)));hold on;
        plot(1:length(marker(three4thIndex,2,:)),squeeze(marker(three4thIndex,2,:))); legend;
    end

    % Nondimensionalizing
    if(0)
        marker(:,1,:) = 1/delta*marker(:,1,:);
        marker(:,2,:) = 1/delta*(marker(:,2,:)-panelLocationY);
        marker(:,8,:) = 1/pinf*marker(:,8,:);
    end

    if(1)
        % Removing the part that is not flexibile
        if(iterE(caseNumber)==0); iterE(caseNumber) = length(marker(1,1,:)); end
        originalMarker = marker; clear marker;
        marker = originalMarker(:,:,iterS(caseNumber):iterE(caseNumber));
    end
    nTime = length(marker(1,1,:));
    t_norm = fluid.ntec*fluid.dt*fluid.Uinf/delta*(1:nTime);
    t_normAllCases{caseNumber} = t_norm;
    clear defAtMiddle def34th def14th;
    defAtMiddle(:) = marker(midPointIndex,2,:);
    def34th(:) = marker(three4thIndex,2,:);
    def14th(:) = marker(one4thIndex,2,:);
    defAtMiddleAllCases{caseNumber} = defAtMiddle;
    def34thAllCases{caseNumber} = def34th;
    def14thAllCases{caseNumber} = def14th;
    tempAverageDefatMiddle = mean(defAtMiddle);
    maxDef = max(defAtMiddle);
    %minDef = min(defAtMiddle);
    nParts = 2;

    %% getting the peaks and valleys of the oscillation
    [peaks,peaksTimeIndex,valleys,valleysTimeIndex,dampingRatio,freqValues] = ...
        plotDeflectionAndDampingRatio(folderNames{caseNumber},imageFormat,defAtMiddle,t_norm,nParts);
    saveas(gcf,"user_dampingRatio_"+folderNames{caseNumber},imageFormat);
    dampingRatioAllCases{caseNumber} = dampingRatio;
    freqValuesAllCases{caseNumber} = freqValues;
    %%
    %=============================================================================
    % Selecting time zone for analysis
    %=============================================================================
    tempIndexS = 6;
    tempIndexE = length(valleysTimeIndex);

    if(tempIndexS <= tempIndexE-5)
        startValue = valleysTimeIndex(tempIndexS);
        endValue = valleysTimeIndex(tempIndexE);
    else
        tempIndexS = max([1,tempIndexE-5]);
        disp("Considering last 5 peaks or all peaks if less than 5")
        startValue = valleysTimeIndex(tempIndexS);
        endValue = valleysTimeIndex(tempIndexE);
    end
    figure(12323);
    plot(t_norm,defAtMiddle,'--r');
    hold on;
    plot(t_norm(startValue:endValue),defAtMiddle(startValue:endValue),'k');
    legend("raw data","zone selected for analysis")
    %%
    %=============================================================================
    % Mid point deflection
    %=============================================================================
    plot(t_norm,defAtMiddle,'-k');
    % text(fluid.dt*peaksTimeIndex,peaks,'*')
    % text(fluid.dt*valleysTimeIndex,valleys,'*')
    titleText = sprintf('T_w/T_r = %s and C_a = %.2e (case:%d)', num2str(round(solid.E/T_r*100)/100),CauchyNumber(caseNumber),caseNumber);
    grid on
    %     title(['(a) Oscillation at the middle of the panel for ',titleText],'fontweight','normal')
    xlabel('$tU/a$','Interpreter','latex'); ylabel('$y/a$','Interpreter','latex')
    set(gca, 'FontName', 'Times')
    set(gcf,'Position',[100,100,850,225]);
    saveas(gcf,"user_OnlyDeflection_"+folderNames{caseNumber},imageFormat)

    %% Plotting time-colored data
    averageDefatMiddle(caseNumber) = mean(defAtMiddle(startValue:endValue));
    nPeaks = length(peaksTimeIndex);
    nValleys = length(valleysTimeIndex);
    peaksTimeIndexAllCases{caseNumber} = peaksTimeIndex;
    valleysTimeIndexAllCases{caseNumber} = valleysTimeIndex;

    J = customcolormap([0 0.5 1], {'#00FF00','#ff0000','#0000FF'},nTime); color = colormap(J);
    figure(7677);clf;
    rowsPlot = 3; colPlot = 2;

    %=============================================================================
    % Mid point deflection with timed colorbar
    %=============================================================================
    tiledplot = tiledlayout(rowsPlot,colPlot,'TileSpacing','Compact','Padding','Compact');
    nexttile([1 2])
    step = 2;
    for i = 1:step:nTime-step
        plot(t_norm(i:i+step),defAtMiddle(i:i+step),'Color',color(i,:));
        hold on
    end
    % text(fluid.dt*peaksTimeIndex,peaks,'*')
    % text(fluid.dt*valleysTimeIndex,valleys,'*')
    disp(['T recovery is ',num2str(T_r)])
    titleText = sprintf('T_w/T_r = %s and C_a = %.2e (case:%d)', num2str(round(fluid.Twall/T_r*100)/100),CauchyNumber(caseNumber),caseNumber);
    grid on
    xlabel('$tU/a$','Interpreter','latex'); ylabel('$y/a$','Interpreter','latex')
    set(gca, 'FontName', 'Times')


    %=============================================================================
    % Panel deflection at Max position
    %=============================================================================
    nexttile;
    for i=1:length(peaks)
        plot(marker(:,1,peaksTimeIndex(i)), marker(:,2,peaksTimeIndex(i)),'Color',color(peaksTimeIndex(i),:))
        hold on
    end
    grid on
    %     title('(b) Panel profiles at min deflection','fontweight','normal')
    xlabel('x/a'); ylabel('y/a')
    axis([xMarkerFlexS/delta xMarkerFlexE/delta -inf inf])
    set(gca, 'FontName', 'Times')

    %=============================================================================
    % Panel deflection at Min position
    %=============================================================================
    nexttile;
    for i=1:length(valleys)
        plot(marker(:,1,valleysTimeIndex(i)),marker(:,2,valleysTimeIndex(i)),'Color',color(valleysTimeIndex(i),:))
        hold on
    end
    grid on
    %     title('(c) Panel profiles at max deflection','fontweight','normal')
    xlabel('x/a');ylabel('y/a')
    axis([xMarkerFlexS/delta xMarkerFlexE/delta -Inf Inf])
    set(gca, 'FontName', 'Times')
    set(gcf, 'Position',  [100, 100, 800, 450])

    %=============================================================================
    % Pressure at min deflection
    %=============================================================================
    nexttile;
    for i=1:length(peaks)
        plot(marker(:,1,peaksTimeIndex(i)), marker(:,8,peaksTimeIndex(i)),'Color',color(peaksTimeIndex(i),:))
        hold on
    end
    grid on
    %     title("(d) Pressure profiles at min deflection",'fontweight','normal')
    xlabel('x/a'); ylabel('P/P_\infty')
    axis([xMarkerS/delta xMarkerE/delta 0.7 1.42])
    set(gca, 'FontName', 'Times')

    %=============================================================================
    % Pressure at max deflection
    %=============================================================================
    nexttile;
    for i=1:length(valleys)
        plot(marker(:,1,valleysTimeIndex(i)),marker(:,8,valleysTimeIndex(i)),'Color',color(valleysTimeIndex(i),:))
        hold on
    end
    %     title("(e) Pressure profiles at max deflection",'fontweight','normal')
    xlabel('x/a'); ylabel('P/P_\infty')
    axis([xMarkerS/delta xMarkerE/delta 0.7 1.42])
    grid on
    set(gca, 'FontName', 'Times')
    set(gcf, 'Position',  [100, 100, 800, 450])

    saveas(tiledplot,"user_n_1_"+num2str(caseNumber)+"_Def_Osci_Pr_"+folderNames{caseNumber},imageFormat)
    %     disp(['SAVED: n_1_',num2str(caseNumber),'_Def_Osci_Pr_',folderNames{caseNumber}],imageFormat)

    %% =============================================================================
    % PSD of Oscillation and Pressure
    %=============================================================================

    %     disp('ACHTUNG!! please check timestep')
    if(solid.thermoelasticFlag==1)
        variableListForFFT = [2 8 9 11];
    else
        variableListForFFT = [2 8 9];
    end

    nColsFFT = length(variableListForFFT);
    nRowsFFT = 1;

    variable = variableListForFFT(1); % 2 for deflection, 8 for pressure
    signalFFT2(:) = marker(midPointIndex,variable,startValue:endValue);
    signalFFT2 = signalFFT2 - mean(signalFFT2);
    t_FFT = t_norm(startValue:endValue);
    [f_2,xt_2] = x_f_PSD_2(t_FFT,signalFFT2); % Doing FFT here - values are nondimensionalized inside
    [MaxMag,IndexTemp] = max(xt_2);
    freqValue(caseNumber) = f_2(IndexTemp);
    figure(1927); clf; plot(f_2,xt_2);
    titleText1 = sprintf('case %d',caseNumber); title(titleText1);
    text(f_2(IndexTemp),MaxMag,['\leftarrow f = ',num2str(f_2(IndexTemp))]);
    saveas(gcf,"user_n_2_Case_"+(caseNumber)+"_PSD_2_"+folderNames{caseNumber},imageFormat);

    [sigr, coeffiecients] = fit_damped_sinewave(signalFFT2);
    coeffiecients(2) = fluid.ntec*fluid.dt*coeffiecients(2)*delta/fluid.Uinf;
    coeffiecients(3) = fluid.ntec*fluid.dt*coeffiecients(3)*delta/fluid.Uinf;
    dampingRatio2 = coeffiecients(2)/sqrt(coeffiecients(2)^2+coeffiecients(3)^2);
    ffidd = fopen(fullfile(dirName{caseNumber},folderNames{caseNumber},"OscillationData.dat"),'W');
    sprintf('A = %f, alpha = %f, w/(2pi) = %f, phi = %f, xi = %f',coeffiecients,dampingRatio2)
    fprintf(ffidd,'A = %f, alpha = %f, w/(2pi) = %f, phi = %f, xi = %f',coeffiecients,dampingRatio2);
    fclose(ffidd);
    fprintf(fidDampedSinosoid,'%s,',folderNames{caseNumber});
    fprintf(fidDampedSinosoid,'%f,%f,%f,%f,%f\n',coeffiecients,dampingRatio2);
    coeffiecientsAllCases(caseNumber,1:4) = coeffiecients;
    coeffiecientsAllCases(caseNumber,5) = dampingRatio2;
    saveas(gcf,"user_fitted_damped_sinewave_"+caseNumber+folderNames{caseNumber},imageFormat);
    clear signalFFT2 f1_2 xt1_2 f_2 xt_2 f_2 xt_2 IndexTemp

    figure(7676); clf;
    tiledplotFFT = tiledlayout(nRowsFFT,nColsFFT,'TileSpacing','Compact','Padding','Compact');
    x = marker(:,1,1);
    for ii = 1:nColsFFT
        nexttile;
        hold on;
        variable = variableListForFFT(ii); % 2 for deflection, 8 for pressure
        signalFFT(:,:) = marker(:,variable,startValue:endValue);
        [f,xt] = x_f_PSD(t_FFT,signalFFT); % Doing FFT here - values are nondimensionalized inside
        titleText1 = sprintf('case %d',caseNumber);
        imagesc('XData',x,'YData',f,'CData',xt');

        colormap(flipud(hot)); colorbar; set(gca,'ColorScale','log')

        if (variable==8)
            %             title(['(b) Pressure variation PSD for ',titleText1],'Interpreter','latex','fontweight','normal')
            %             caxis([0 10])
            %             ylim([0 0.01])
            %             xlim([xMarkerS/delta xMarkerE/delta])
        elseif(variable==2)
            %             title(['(a) Panel oscillation PSD for ',titleText1],'Interpreter','latex','fontweight','normal')
            %             caxis([0 80])
            %             ylim([0 0.01])
            %             xlim([xMarkerS/delta xMarkerE/delta])
        elseif(variable==10)
            %             title(['Shear stress variation PSD for ',titleText1],'Interpreter','latex','fontweight','normal')
            %caxis([0 5])
            %ylim([0 5])
            %             xlim([xMarkerS/delta xMarkerE/delta])
        end
        hold on
        xline(min(xlim))
        xline(max(xlim))
        yline(max(ylim))
        ylabel('$fa/U$','Interpreter','latex')
        xlabel('$x/a$','Interpreter','latex')
        set(gca, 'FontName', 'Times')
        grid on;        grid minor
    end
    set(gcf, 'Position',  [100, 100, 800, 500])
    saveas(tiledplotFFT,"user_n_2_Case_"+caseNumber+"_PSD_"+folderNames{caseNumber},imageFormat)
    %     disp(['SAVED: n_2_Case_',num2str(caseNumber),'_PSD_',folderNames{caseNumber}],imageFormat)
    clear signalFFT t_FFT

    %% =============================================================================
    % Avg on Valley and Peaks vs Instantaneous
    %=============================================================================

    tempIndexS = 6;
    tempIndexE = length(valleysTimeIndex);

    if(tempIndexS <= tempIndexE-5)
        startValue = valleysTimeIndex(tempIndexS);
        endValue = valleysTimeIndex(tempIndexE);
    else
        tempIndexS = max([1,tempIndexE-5]);
        disp("Considering last 5 peaks or all peaks if less than 5")
        startValue = valleysTimeIndex(tempIndexS);
        endValue = valleysTimeIndex(tempIndexE);
    end

    figure(12323); clf;
    plot(t_norm,def34th/0.002,':r');
    hold on;
    plot(t_norm(startValue:endValue),def34th(startValue:endValue)/0.002,'k');
    legend("deflection","averaging window")
    set(gcf, 'Position',  [100, 100, 300, 250])
    ylabel('$y/h$','Interpreter','latex')
    xlabel('$tU/a$','Interpreter','latex')
    set(gca, 'FontName', 'Times', 'ycolor','k', 'FontSize',10)

    %% Plot averaging data
    variableForAvg = 8;
    figure(10089); clf;
    hold on;
    yyaxis left;

    clear avgMarker
    avgMarker=0;
    count=0;
    % average over 1st valleys to last valleys
    for i=startValue:endValue
        avgMarker = avgMarker + marker(:,:,i);
        count=count+1;
    end
    avgMarker = avgMarker/count;
    plot(marker(:,1,1),avgMarker(:,variableForAvg));
    if(solid.thermoelasticFlag~=1);  avgMarker(:,11) = 0; end
    multipleCasesPr(:,:,caseNumber) = avgMarker(:,:);
    avgMaxMarker=0;
    count=0;

    for i=1:length(valleys)% average only at the valleys
        avgMaxMarker = avgMaxMarker + marker(:,:,valleysTimeIndex(i));
        count=count+1;
    end
    avgMaxMarker = avgMaxMarker/count;
    plot(marker(:,1,1),avgMaxMarker(:,variableForAvg));

    avgMinMarker=0;
    count=0;
    for i=1:length(peaks)% average only at the peaks
        avgMinMarker = avgMinMarker + marker(:,:,peaksTimeIndex(i));
        count=count+1;
    end
    avgMinMarker = avgMinMarker/count;
    plot(marker(:,1,1),avgMinMarker(:,variableForAvg));

    %     title(['Pressure profile for ',titleText],'fontweight','normal')
    ylabel('$P/P_\infty$','Interpreter','latex')
    xlabel('$x/a$','Interpreter','latex')
    axis([xMarkerS/delta xMarkerE/delta 0.7 1.1])

    [~,sDefCurve] = min(abs(xx-xMarkerFlexS));
    [~,eDefCurve] = min(abs(xx-xMarkerFlexE));

    % showing deflected part
    yyaxis right;
    plot(avgMinMarker(sDefCurve:eDefCurve,1),avgMinMarker(sDefCurve:eDefCurve,2),avgMaxMarker(sDefCurve:eDefCurve,1),avgMaxMarker(sDefCurve:eDefCurve,2))
    ylabel('$y/a$','Interpreter','latex')
    % axis([xMarkerS/delta xMarkerE/delta 0.5/delta 10/delta])
    xline(xMarkerFlexS/delta,'--','Color',[0.4 0.4 0.4]);
    xline(xMarkerFlexE/delta,'--','Color',[0.4 0.4 0.4]);
    grid on
    legend('Flexible(Avg)','Flexible(Min)','Flexible(Max)','Right Axis: Min Def','Right Axis: Max Def','Location','northwest')
    set(gca, 'FontName', 'Times')
    set(gcf, 'Position',  [100, 100, 450, 250])
    saveas(gca,"user_n_5_"+(caseNumber)+"_PressureProfile_"+folderNames{caseNumber},imageFormat)
    disp(['SAVED: user_n_5_',num2str(caseNumber),'_PressureProfile_',folderNames{caseNumber}])

    %% =============================================================================
    % Avg profiles vs Instanteneous
    %=============================================================================

    figure(834998); clf;
    tiledplotAvg = tiledlayout(2,1,'TileSpacing','Compact','Padding','Compact');
    nexttile;
    plot(t_norm,defAtMiddle,'-k');hold on;
    plot(t_norm(valleysTimeIndex(tempIndexS+1)),defAtMiddle(valleysTimeIndex(tempIndexS+1)),'*r');
    plot(t_norm(peaksTimeIndex(tempIndexS)),defAtMiddle(peaksTimeIndex(tempIndexS)),'*b');
    nexttile;
    yyaxis left;
    hold on;
    plot(marker(:,1,1),avgMarker(:,variableForAvg),'-k');
    plot(marker(:,1,1),marker(:,variableForAvg,valleysTimeIndex(tempIndexS+1)),'-r');
    plot(marker(:,1,1),marker(:,variableForAvg,peaksTimeIndex(tempIndexS+1)),'-b');

    %     title(['Pressure profile for ',titleText],'fontweight','normal')
    ylabel('$P/P_\infty$','Interpreter','latex')
    xlabel('$x/a$','Interpreter','latex')
    axis([xMarkerS/delta xMarkerE/delta 0.7 1.1])
    set(gca, 'FontName', 'Times', 'ycolor','k')

    % showing deflected part
    yyaxis right;
    [~,sDefCurve] = min(abs(xx-xMarkerFlexS));
    [~,eDefCurve] = min(abs(xx-xMarkerFlexE));
    plot(avgMinMarker(sDefCurve:eDefCurve,1),avgMinMarker(sDefCurve:eDefCurve,2),'-.b',avgMaxMarker(sDefCurve:eDefCurve,1),avgMaxMarker(sDefCurve:eDefCurve,2),'-.r')
    ylabel('$y/a$','Interpreter','latex')
    % axis([xMarkerS/delta xMarkerE/delta 0.5/delta 2.5/delta])
    xline(xMarkerFlexS/delta,'--','Color',[0.4 0.4 0.4]);
    xline(xMarkerFlexE/delta,'--','Color',[0.4 0.4 0.4]);
    legend('Flexible(Avg)','Flexible(Max)','Flexible(Min)','Panel(Min)','Panel(Max)','Location','northwest')
    grid on
    set(gca, 'FontName', 'Times', 'ycolor','k')
    set(gcf, 'Position',  [100, 100, 400, 250])
    saveas(gca,"user_n_5_"+(caseNumber)+"_PressureProfile_AvgVsInst"+folderNames{caseNumber},imageFormat)
    %     disp(['SAVED: n_5_',num2str(caseNumber),'_PressureProfile_AvgVsInst',folderNames{caseNumber}],imageFormat)

end
fclose(fidDampedSinosoid);
%% Saving Damped Sinosoid information
save('freqValue_user','freqValue');
save('coeffiecientsAllCases_user','coeffiecientsAllCases','folderNames','dirName');
save('multipleCasesPr_user','multipleCasesPr','folderNames','dirName');
save('defAtMiddleAllCases_user','defAtMiddleAllCases','t_normAllCases','folderNames','dirName');
save('def34thAllCases_user','def34thAllCases','t_normAllCases','folderNames','dirName');
save('def14thAllCases_user','def14thAllCases','t_normAllCases','folderNames','dirName');
save('valleysTimeIndexAllCases_user','valleysTimeIndexAllCases','folderNames','dirName')

%%
load('coeffiecientsAllCases_user','coeffiecientsAllCases','folderNames','dirName');

nVar = length(coeffiecientsAllCases(1,:));
coefNames{1} = 'amplitude';
coefNames{2} = 'decayRate';
coefNames{3} = 'frequency';
coefNames{4} = 'phase';
coefNames{5} = 'dampingRatio';

coefNamesX{1} = '$A/a$';
coefNamesX{2} = '$\lambda a /U$';
coefNamesX{3} = '$fa /U$';
coefNamesX{4} = '$\phi$';
coefNamesX{5} = '$\zeta$';

titleName{1} = '(a)';
titleName{2} = '(b)';
titleName{3} = '(c)';
titleName{4} = '(d)';
titleName{5} = '(d)';

colorS = ['r';'b';'g'];

figure;
clear xTemp yTemp

nCases = length(folderNames);

for iVar = 1:nVar
    figure; clf;
    plot( 1:nCases,coeffiecientsAllCases(:,iVar));

    Legend = cell(nCases,1);
    for ic = 1:nCases
        Legend{ic} = ['C_a = ', sprintf('%.2e',CauchyNumber(ic))];
    end

    grid on
    %xlim([min(n_T)/T_r max(n_T)/T_r])
    xlabel('Case number','Interpreter','latex')
    ylabel(coefNamesX{iVar},'Interpreter','latex')
    title(titleName{iVar},'Interpreter','latex')
    set(gca, 'FontName', 'Times')
    set(gcf, 'Position',  [100, 100, 175, 250])
    imageName = ['user_dampedSinosoid',coefNames{iVar}];
    if(iVar~=5)
        saveas(gcf,imageName,imageFormat)
    end
end
set(gcf, 'Position',  [100, 100, 315, 250])
legend(Legend,'location','eastoutside')
saveas(gcf,imageName,imageFormat)

%% Multiple Cases - Def, Pressure, Shear

load('defAtMiddleAllCases_user','defAtMiddleAllCases','t_normAllCases','folderNames','dirName');
load('valleysTimeIndexAllCases_user','valleysTimeIndexAllCases','folderNames','dirName');
load('multipleCasesPr_user','multipleCasesPr','folderNames','dirName');

casesTocomp = [1 6];
nCases = length(casesTocomp);
figure(850);
variableList = [2 8 10]; % Pressure
color = ['r';'g';'b'];

Legend = cell(nCases,1);
for iT = 1:nCases
    Legend{iT} = ['Case ', num2str(iT)];
end

figure(89763);clf;
figure(89764);clf;
figure(89765);clf;

% plotting timeseries data - showing the avging window
figure(89763);
for i = 1:nCases
    caseNumber = casesTocomp(i);
    hold on
    plot(t_normAllCases{caseNumber},(defAtMiddleAllCases{caseNumber}-1)/delta,'--','Color',color(i));
end

for i = 1:nCases
    caseNumber = casesTocomp(i);
    valleysTimeIndex = valleysTimeIndexAllCases{caseNumber};
    tempIndexS = 6;
    tempIndexE = 12;

    if(tempIndexE<= length(valleysTimeIndex))
        startValue = valleysTimeIndex(tempIndexS);
        endValue = valleysTimeIndex(tempIndexE);
    else
        tempIndexS = length(valleysTimeIndex) - 5;
        disp("Considering last 5 peaks")
        tempIndexE = length(valleysTimeIndex);
        startValue = valleysTimeIndex(tempIndexS);
        endValue = valleysTimeIndex(tempIndexE);
    end
    hold on
    plot(t_normAllCases{caseNumber}(startValue:endValue),(defAtMiddleAllCases{caseNumber}(startValue:endValue)-1)/delta,'Color',color(i));
end
ylabel('$y/a$','Interpreter','latex')
%     title(['(a) Panel oscillation'],'Interpreter','latex','fontweight','normal');
xlabel('$tU/a$','Interpreter','latex')
set(gca, 'FontName', 'Times')
Legend = cell(nCases,1);
for iT = 1:nCases
    Legend{iT} = ['Case ', num2str(nCases)];
end
legend(Legend,'Location','eastoutside')
set(gcf, 'Position',  [100, 100, 600, 250])
box on; grid on;
imageName = ['user_multipleCases_E_Osci',num2str(caseNumber),'_Variable_',num2str(0)];
saveas(gcf,imageName,imageFormat);

% Avg Pressure and cf
for i = 1:nCases
    caseNumber = casesTocomp(i);
    [fluid] = loadCaseData(dirName{caseNumber},folderNames{caseNumber});
    dynamicPressure = 0.5*rhoinf*fluid.Uinf^2;
    mu = rhoinf*fluid.Uinf*1/fluid.Re;
    figure(89764); % plotting pressure
    hold on
    avgMarker(:,:) = multipleCasesPr(:,:,caseNumber);
    plot(marker(:,1,1)/delta,avgMarker(:,8),'Color',color(i));
    ylabel('$P/P_\infty$','Interpreter','latex')
    ylim([0.8 1.9])
    xlabel('$x/a$','Interpreter','latex')
    set(gca, 'FontName', 'Times','ycolor','b')

    figure(89765); % plotting cf
    hold on
    plot(marker(:,1,1)/delta,mu*avgMarker(:,10)/dynamicPressure,'Color',color(i));
    yline(0);
    ylabel('$\tau/\rho_\infty U^2$','interpreter','latex')
    xlabel('$x/a$','Interpreter','latex')
    set(gca, 'FontName', 'Times','ycolor','b')
    %         ylim([0.8 1.9])
end
% xlim([xMarkerS/delta xMarkerE/delta])

clear avgDiff

figure(89764);
yyaxis right;
avgDiff(:) = abs(multipleCasesPr(:,8,1)-multipleCasesPr(:,8,end));
plot(marker(:,1,1)/delta,avgDiff(:),'--k');
ylabel('$a P/ P_\infty$','Interpreter','latex')
%ylim([0 0.3*4])
grid on
%xlim([xMarkerS/delta xMarkerE/delta])
set(gca, 'FontName', 'Times','ycolor','k')
set(gcf, 'Position',  [500, 100, 400, 250])
legend(Legend,'Location','best')
imageName = ['user_multipleCases_E_Pressure',num2str(caseNumber),'_Variable_',num2str(0)];
box on; grid on;
saveas(gcf,imageName,imageFormat);

figure(89765);
yyaxis right;
avgDiff(:) = abs(multipleCasesPr(:,10,1)-multipleCasesPr(:,10,end));
plot(marker(:,1,1)/delta,avgDiff(:),'--k');
ylabel('$a\tau / \rho_\infty U^2$','Interpreter','latex')
%ylim([0 40])
grid on
set(gca, 'FontName', 'Times','ycolor','k')
set(gcf, 'Position',  [900, 100, 400, 250])
legend(Legend,'Location','best')
imageName = ['user_multipleCases_E_Tau',num2str(0),'_Variable_',num2str(0)];
box on; grid on;
saveas(gcf,imageName,imageFormat);

%% Writing
fileID = fopen('peaks_user.txt','wt');
for caseNumber = 1:totalNumberOfCases
    nPeaks = length(peaksTimeIndexAllCases);
    fprintf(fileID,'%s %d \n',folderNames{caseNumber},nPeaks);
    fprintf(fileID,'%d ',iterS(caseNumber)+peaksTimeIndexAllCases{caseNumber}-1);
    fprintf(fileID,'\n');
end
fclose(fileID);

fileID = fopen('valleys_user.txt','wt');
for caseNumber = 1:totalNumberOfCases
    nValleys = length(valleysTimeIndexAllCases{caseNumber});
    fprintf(fileID,'%s %d \n',folderNames{caseNumber},nValleys);
    fprintf(fileID,'%d ',iterS(caseNumber)+valleysTimeIndexAllCases{caseNumber}-1);
    fprintf(fileID,'\n');
end
fclose(fileID);
%
fileID = fopen('peaks2_user.txt','wt');
for caseNumber = 1:totalNumberOfCases
    nPeaks = length(peaksTimeIndexAllCases{caseNumber});
    fprintf(fileID,'%d ',peaksTimeIndexAllCases{caseNumber});
    fprintf(fileID,'\n');
end
fclose(fileID);

fileID = fopen('valleys2_user.txt','wt');
for caseNumber = 1:totalNumberOfCases
    nValleys = length(valleysTimeIndexAllCases{caseNumber});
    fprintf(fileID,'%d ',valleysTimeIndexAllCases{caseNumber});
    fprintf(fileID,'\n');
end
fclose(fileID);
