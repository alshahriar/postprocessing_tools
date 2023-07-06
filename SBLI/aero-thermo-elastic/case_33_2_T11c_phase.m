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
    elseif(1)
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
averageDefatMiddle = zeros(1,totalNumberOfCases);
CauchyNumber = zeros(1,totalNumberOfCases);
freqValue = zeros(1,totalNumberOfCases);



%%
for caseNumber = 2:2
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
    if(0)
        figure;
        % quick check on the time series of deflection
        plot(1:length(marker(midPointIndex,2,:)),squeeze(marker(midPointIndex,2,:)));hold on;
        plot(1:length(marker(three4thIndex,2,:)),squeeze(marker(three4thIndex,2,:))); legend;
    end

    % Nondimensionalizing
    if(1)
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
        plotDeflectionAndDampingRatio(folderNames{caseNumber},imageFormat,def34th,t_norm,nParts);
    saveas(gcf,"user_dampingRatio_"+folderNames{caseNumber},imageFormat);
    dampingRatioAllCases{caseNumber} = dampingRatio;
    freqValuesAllCases{caseNumber} = freqValues;
    %% Doing phase averaging
    %=============================================================================
    % Selecting time zone for analysis
    %=============================================================================

    timeS = 169.;
    timeE = 369;
    [~,startValue] = min(abs(t_norm-timeS));
    [~,endValue] = min(abs(t_norm-timeE));

    %mrktrim = marker(:,:,startValue:endValue);

    figure(12323);
    plot(t_norm,defAtMiddle,'--r');
    hold on;
    plot(t_norm(startValue:endValue),defAtMiddle(startValue:endValue),'k');
    legend("raw data","zone selected for analysis")


    nValleys_phase = 24;
    disp("valleys for phase average"+nValleys_phase);
    phaseGap = round((endValue-startValue)/nValleys_phase);
    marker_phase = zeros(length(marker(:,1,1)),length(marker(1,:,1)),phaseGap);
    endValue = startValue + nValleys_phase*phaseGap;
    for i = 1:phaseGap
        js = startValue+i-1 + phaseGap*((1:nValleys_phase)-1);
        for j = js
            %disp(j)
            marker_phase(:,:,i) = marker_phase(:,:,i) + marker(:,:,j);
        end
    end
    marker_phase = marker_phase/nValleys_phase;

    figure;
    for i=1:phaseGap
        clf;
        plot(marker_phase(:,1,i), marker_phase(:,2,i))
        ylim([-0.02 0.02])
        drawnow;
        pause(0.05);
    end
    %     figure; hold on;
    %     for i=1:phaseGap
    %         clf;
    %         plot(marker_phase(:,1,i), marker_phase(:,9,i))
    %         %ylim([-0.02 0.02])
    %         drawnow;
    %         pause(0.05);
    %     end

    t_ph(:,1) = t_norm(1:phaseGap);
    d_mid_ph(:,1) = marker_phase(midPointIndex,2,:);
    d_34_ph(:,1) = marker_phase(three4thIndex,2,:);
    %%
    %=============================================================================
    % Mid point deflection
    %=============================================================================
    figure;
    plot(t_ph,d_mid_ph,'-k');
    titleText = sprintf('T_w/T_r = %s and C_a = %.2e (case:%d)', num2str(round(solid.E/T_r*100)/100),CauchyNumber(caseNumber),caseNumber);
    grid on
    %     title(['(a) Oscillation at the middle of the panel for ',titleText],'fontweight','normal')
    xlabel('$tU/a$','Interpreter','latex'); ylabel('$y/a$','Interpreter','latex')
    set(gca, 'FontName', 'Times')
    set(gcf,'Position',[100,100,850,225]);
    saveas(gcf,"user_OnlyDeflection_"+folderNames{caseNumber},imageFormat)

    %% Plotting time-colored data
    averageDefatMiddle(caseNumber) = mean(d_mid_ph);
    nLines = 5;
    snapshot_time_index = round(linspace(1,phaseGap,nLines));

    J = customcolormap([0 0.5 1], {'#00FF00','#FF0000','#0000FF'},phaseGap); color = colormap(J);
    figure(7677);clf;
    rowsPlot = 3; colPlot = 2;
    %=============================================================================
    % Mid point deflection with timed colorbar
    %=============================================================================
    tiledplot = tiledlayout(rowsPlot,colPlot,'TileSpacing','Compact','Padding','Compact');
    nexttile([1 2])
    step = 2;
    for i = 1:step:phaseGap-step
        plot(t_ph(i:i+step),d_mid_ph(i:i+step),'Color',color(i,:));
        hold on
    end

    for i = 1:step:phaseGap-step
        plot(t_ph(i:i+step),d_34_ph(i:i+step),'--','Color',color(i,:));
        hold on
    end
    disp(['T recovery is ',num2str(T_r)])
    titleText = sprintf('T_w/T_r = %s and C_a = %.2e (case:%d)', num2str(round(fluid.Twall/T_r*100)/100),CauchyNumber(caseNumber),caseNumber);
    grid on
    xlabel('$tU/a$','Interpreter','latex'); ylabel('$y/a$','Interpreter','latex')
    set(gca, 'FontName', 'Times')

    %=============================================================================
    % Panel deflection at Max position
    %=============================================================================
    nexttile;
    for i=1:nLines
        plot(marker_phase(:,1,snapshot_time_index(i)), marker_phase(:,2,snapshot_time_index(i)),'Color',color(snapshot_time_index(i),:))
        hold on
    end
    grid on
    %     title('(b) Panel profiles at min deflection','fontweight','normal')
    xlabel('$x/a$','Interpreter','latex'); ylabel('$y/a$','Interpreter','latex')
    xlim([xMarkerFlexS/delta xMarkerFlexE/delta])
    %axis([xMarkerFlexS/delta xMarkerFlexE/delta -inf inf])
    set(gca, 'FontName', 'Times')

    %=============================================================================
    % Pressure
    %=============================================================================
    nexttile;
    for i=1:nLines
        plot(marker_phase(:,1,snapshot_time_index(i)),marker_phase(:,8,snapshot_time_index(i)),'Color',color(snapshot_time_index(i),:))
        hold on
    end
    grid on
    %     title('(c) Panel profiles at max deflection','fontweight','normal')
    xlabel('$x/a$','Interpreter','latex');ylabel('$P/P_\infty$','Interpreter','latex')
    xlim([xMarkerFlexS/delta xMarkerFlexE/delta])
    %axis([xMarkerFlexS/delta xMarkerFlexE/delta -Inf Inf])
    set(gca, 'FontName', 'Times')
    set(gcf, 'Position',  [100, 100, 900, 450])

    %=============================================================================
    % v velocity
    %=============================================================================
    nexttile;
    for i=1:nLines
        plot(marker_phase(:,1,snapshot_time_index(i)), marker_phase(:,6,snapshot_time_index(i)),'Color',color(snapshot_time_index(i),:))
        hold on
    end
    grid on
    %     title("(d) Pressure profiles at min deflection",'fontweight','normal')
    xlabel('$x/a$','Interpreter','latex'); ylabel('$v/c_\infty$','Interpreter','latex')
    xlim([xMarkerFlexS/delta xMarkerFlexE/delta])
    %axis([xMarkerFlexS/delta xMarkerFlexE/delta 0.7 1.42])
    set(gca, 'FontName', 'Times')

    %=============================================================================
    % Heatflux
    %=============================================================================
    nexttile;
    for i=1:nLines
        plot(marker_phase(:,1,snapshot_time_index(i)),marker_phase(:,11,snapshot_time_index(i)),'Color',color(snapshot_time_index(i),:))
        hold on
    end
    %     title("(e) Pressure profiles at max deflection",'fontweight','normal')
    %xlabel('x/a'); ylabel('Nu')
    xlabel('$x/a$','Interpreter','latex'); ylabel('$Nu$','Interpreter','latex')
    %axis([xMarkerFlexS/delta xMarkerFlexE/delta 0.7 1.42])
    xlim([xMarkerFlexS/delta xMarkerFlexE/delta])
    grid on
    set(gca, 'FontName', 'Times')
    set(gcf, 'Position',  [100, 100, 800, 450])

    saveas(tiledplot,"phase_avg_"+num2str(caseNumber)+"_Def_Osci_Pr_"+folderNames{caseNumber},'fig')
    %     disp(['SAVED: n_1_',num2str(caseNumber),'_Def_Osci_Pr_',folderNames{caseNumber}],imageFormat)

    %% xt plot

    % xt plot of variables
    variableList_xt = [2 4 6 8 9 11];
    variableNames_xt = 'yrvPTQ';

    colormapType = [1 2 1 2 2 2];

    J = customcolormap([0 0.5 1], {'#0000FF','#FFFFFF','#FF0000'},100);
    BlWtGr = colormap(J);

    nVar_xt = length(variableList_xt);
    xxt(:) = marker_phase(:,1,1);
    [~,indexFlexS] = min(abs(xxt-xMarkerFlexS));
    [~,indexFlexE] = min(abs(xxt-xMarkerFlexE));
    %[~,indexFlexS] = min(abs(xxt+14));
    %[~,indexFlexE] = min(abs(xxt-14));

    clear marker_n xtV
    marker_n(:,:,:) = marker_phase(indexFlexS+2:indexFlexE-2,:,:);
   
    if(solid.thermoelasticFlag==1)
        xtV(:,:) = marker_n(:,11,:);
        %minForQ = min(min(xtV));
        %maxForQ = max(max(xtV));
        minForQ = -1e-3;
        maxForQ = 1e-3;
    else
        minForQ = 1;
        maxForQ = 2;
    end

    xtV(:,:) = marker_n(:,9,:)-0;
    minForT = min(min(xtV));
    maxForT = max(max(xtV));
    clim1 = [-0.01 0.6 -5e-3 1.0 minForT minForQ];
    clim2 = [ 0.01 1.2  5e-3 1.4 maxForT maxForQ];
    if(solid.thermoelasticFlag==1)
        nVar_xt_temp = nVar_xt;
    else
        nVar_xt_temp =  nVar_xt - 1;
    end

    for iVar = 1:nVar_xt_temp
        figure;clf;
        varValue = variableList_xt(iVar);
        clear xx xtV;
        xtV(:,:) = marker_n(:,varValue,:);
        xx(:) = marker_n(:,1,1);
        [xxx,tt] = meshgrid(xx,t_ph);
        contourf(xxx,tt,xtV',500,"edgeColor","none")
        if(colormapType(iVar)==1)
            colormap(BlWtGr);
        else
            colormap("jet");
        end
        colorbar;
        set(gca, 'FontName', 'Times', 'ycolor','k')
        set(gcf, 'Position',  [100, 500, 200, 250])

        %ylim([0 300])
        caxis([clim1(iVar) clim2(iVar)])
        xlim([xMarkerFlexS xMarkerFlexE]/delta)
        ylabel('$tU/a$','Interpreter','latex'); xlabel('$x/a$','Interpreter','latex')
        saveas(gca,"phase_avg_user_xt_"+variableNames_xt(iVar)+"_"+folderNames{caseNumber},'fig')
        %saveas(gca,"user_xt_"+variableNames_xt(iVar)+"_"+folderNames{caseNumber},'fig')
    end



end
fclose(fidDampedSinosoid);
