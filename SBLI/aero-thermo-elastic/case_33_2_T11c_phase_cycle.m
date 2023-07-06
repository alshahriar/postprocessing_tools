%% Section 1
goAhead = input('Are you sure?');
clear all

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

    timeS = 286.7;
    timeE = 303.3;
    [~,startValue] = min(abs(t_norm-timeS));
    [~,endValue] = min(abs(t_norm-timeE));

    %mrktrim = marker(:,:,startValue:endValue);

    figure(12323);
    plot(t_norm,defAtMiddle,'--r');
    hold on;
    plot(t_norm(startValue:endValue),defAtMiddle(startValue:endValue),'k');
    legend("raw data","zone selected for analysis")


    cycleGap = endValue - startValue + 1;
    marker_cycle = zeros(length(marker(:,1,1)),length(marker(1,:,1)),cycleGap);
    for i = 1:cycleGap        
        marker_cycle(:,:,i) = marker(:,:,startValue+i-1);
    end

    J = customcolormap([0 0.5 1], {'#00FF00','#FF0000','#0000FF'},cycleGap); color = colormap(J);
    figure;hold on;
    for i = 1:2:cycleGap
        plot(marker_cycle(:,1,i), marker_cycle(:,2,i),'Color',color(i,:));
        hold on
        ylim([-0.02 0.02])
        drawnow;
        pause(0.05);
    end
    clear t_cy d_mid_cy d_34_cy
    t_cy(:,1) = t_norm(startValue:endValue);
    d_mid_cy(:,1) = marker_cycle(midPointIndex,2,:);
    d_34_cy(:,1) = marker_cycle(three4thIndex,2,:);
    %%
    %=============================================================================
    % Mid point deflection
    %=============================================================================
    figure;
    plot(1:cycleGap,d_mid_cy,'-k'); hold on;
    plot(1:cycleGap,d_34_cy,'-r');
    titleText = sprintf('T_w/T_r = %s and C_a = %.2e (case:%d)', num2str(round(solid.E/T_r*100)/100),CauchyNumber(caseNumber),caseNumber);
    grid on
    %     title(['(a) Oscillation at the middle of the panel for ',titleText],'fontweight','normal')
    xlabel('$tU/a$','Interpreter','latex'); ylabel('$y/a$','Interpreter','latex')
    set(gca, 'FontName', 'Times')
    set(gcf,'Position',[100,100,850,225]);
    saveas(gcf,"user_OnlyDeflection_"+folderNames{caseNumber},imageFormat)

    %% Plotting time-colored data
    averageDefatMiddle(caseNumber) = mean(d_mid_cy);
    nLines = 5;
    snapshot_time_index = round(linspace(1,cycleGap,nLines));
    snapshot_time_index = [1 17 41 64 cycleGap];
    nLines = length(snapshot_time_index);
    J = customcolormap([0 0.5 1], {'#00FF00','#FF0000','#0000FF'},cycleGap); color = colormap(J);
    figure(7677);clf;
    rowsPlot = 3; colPlot = 2;
    ttltxt{6} = "(f)";
    ttltxt{7} = "(g)";
    ttltxt{8} = "(h)";
    ttltxt{9} = "(i)";
    ttltxt{10} = "(j)";
    ttltxt{11} = "(k)";    
    %=============================================================================
    % Mid point deflection with timed colorbar
    %=============================================================================
    tiledplot = tiledlayout(rowsPlot,colPlot,'TileSpacing',"tight",'Padding','tight');
    nexttile;
    step = 2;
    for i = 1:step:cycleGap-step
        plot(t_cy(i:i+step),d_mid_cy(i:i+step),'Color',color(i,:));
        hold on
    end

    for i = 1:step:cycleGap-step
        plot(t_cy(i:i+step),d_34_cy(i:i+step),':','Color',color(i,:));
        hold on
    end
    for i=1:nLines
        xline(t_cy(snapshot_time_index(i)),"Color",color(snapshot_time_index(i),:))
    end
    disp(['T recovery is ',num2str(T_r)])
    %titleText = sprintf('T_w/T_r = %s and C_a = %.2e (case:%d)', num2str(round(fluid.Twall/T_r*100)/100),CauchyNumber(caseNumber),caseNumber);
    title(ttltxt{6},'Interpreter','latex')
    xlabel('$tU/a$','Interpreter','latex'); ylabel('$y/a$','Interpreter','latex')
        set(gca, 'FontName', 'Times', 'ycolor','k',"FontSize",9)

    %=============================================================================
    % Panel deflection at Max position
    %=============================================================================
    nexttile;
    for i=1:nLines
        plot(marker_cycle(:,1,snapshot_time_index(i)), marker_cycle(:,2,snapshot_time_index(i)),'Color',color(snapshot_time_index(i),:))
        hold on
    end
    grid on
    %     title('(b) Panel profiles at min deflection','fontweight','normal')
    xlabel('$x/a$','Interpreter','latex'); ylabel('$y/a$','Interpreter','latex')
    xlim([xMarkerFlexS/delta xMarkerFlexE/delta])
    %axis([xMarkerFlexS/delta xMarkerFlexE/delta -inf inf])
    title(ttltxt{7},'Interpreter','latex')
        set(gca, 'FontName', 'Times', 'ycolor','k',"FontSize",9)




    %=============================================================================
    % v velocity
    %=============================================================================
    nexttile;
    for i=1:nLines
        plot(marker_cycle(:,1,snapshot_time_index(i)), marker_cycle(:,6,snapshot_time_index(i)),'Color',color(snapshot_time_index(i),:))
        hold on
    end
    grid on
    %     title("(d) Pressure profiles at min deflection",'fontweight','normal')
    xlabel('$x/a$','Interpreter','latex'); ylabel('$v/c_\infty$','Interpreter','latex')
    xlim([xMarkerFlexS/delta xMarkerFlexE/delta])
    %axis([xMarkerFlexS/delta xMarkerFlexE/delta 0.7 1.42])
    title(ttltxt{9},'Interpreter','latex')
        set(gca, 'FontName', 'Times', 'ycolor','k',"FontSize",9)

    %=============================================================================
    % Pressure
    %=============================================================================
    nexttile;
    avgP = mean(marker_cycle(:,8,:),3);
    for i=1:nLines
        plot(marker_cycle(:,1,snapshot_time_index(i)),marker_cycle(:,8,snapshot_time_index(i))-avgP,'Color',color(snapshot_time_index(i),:))
        hold on
    end
    grid on
    %     title('(c) Panel profiles at max deflection','fontweight','normal')
    xlabel('$x/a$','Interpreter','latex');ylabel('$(P-\overline{P})/P_\infty$','Interpreter','latex')
    xlim([xMarkerFlexS/delta xMarkerFlexE/delta])
    %axis([xMarkerFlexS/delta xMarkerFlexE/delta -Inf Inf])
    title(ttltxt{10},'Interpreter','latex')
        set(gca, 'FontName', 'Times', 'ycolor','k',"FontSize",9)
    set(gcf, 'Position',  [100, 100, 900, 450])
    

    %=============================================================================
    % Density
    %=============================================================================
    nexttile;
    for i=1:nLines
    plot(marker_cycle(:,1,snapshot_time_index(i)), smooth(marker_cycle(:,4,snapshot_time_index(i)),0.07,'loess'),'Color',color(snapshot_time_index(i),:))
        hold on
    end
    grid on
    %     title('(b) Panel profiles at min deflection','fontweight','normal')
    xlabel('$x/a$','Interpreter','latex'); ylabel('$\rho/\rho_\infty$','Interpreter','latex')
    xlim([xMarkerFlexS/delta xMarkerFlexE/delta])
    %axis([xMarkerFlexS/delta xMarkerFlexE/delta -inf inf])
    title(ttltxt{8},'Interpreter','latex')
    set(gca, 'FontName', 'Times', 'ycolor','k',"FontSize",9)



    %=============================================================================
    % Heatflux
    %=============================================================================
    nexttile;
    for i=1:nLines
        plot(marker_cycle(:,1,snapshot_time_index(i)),marker_cycle(:,11,snapshot_time_index(i)),'Color',color(snapshot_time_index(i),:))
        hold on
    end
    %     title("(e) Pressure profiles at max deflection",'fontweight','normal')
    %xlabel('x/a'); ylabel('Nu')
    xlabel('$x/a$','Interpreter','latex'); ylabel('$Nu$','Interpreter','latex')
    %axis([xMarkerFlexS/delta xMarkerFlexE/delta 0.7 1.42])
    xlim([xMarkerFlexS/delta xMarkerFlexE/delta])
    ylim([-0.000 Inf])
    grid on
    title(ttltxt{11},'Interpreter','latex')
            set(gca, 'FontName', 'Times', 'ycolor','k',"FontSize",9)

    set(gcf, 'Position',  [100, 100, 800, 450])

    saveas(tiledplot,"cycle_avg_"+num2str(caseNumber)+"_Def_Osci_Pr_"+folderNames{caseNumber},'fig')
    %     disp(['SAVED: n_1_',num2str(caseNumber),'_Def_Osci_Pr_',folderNames{caseNumber}],imageFormat)

    %% xt plot

    % xt plot of variables
    variableList_xt = [2 6 8 4 11];
    ttltxt{1} = "(a)";
    ttltxt{2} = "(b)";
    ttltxt{3} = "(c)";
    ttltxt{4} = "(d)";
    ttltxt{5} = "(e)";
    variableNames_xt = 'yrvPQ';

    colormapType = [1 2 1 2 2 2];
    J = customcolormap([0 0.5 1], {'#0000FF','#FFFFFF','#FF0000'},100);
    BlWtRd = colormap(J);
    J = customcolormap([0 0.5 1], {'#0000FF','#00FF00','#FF0000'},100);
    BlGrRd = colormap(J);

    nVar_xt = length(variableList_xt);
    xxt(:) = marker_cycle(:,1,1);
    [~,indexFlexS] = min(abs(xxt-xMarkerFlexS));
    [~,indexFlexE] = min(abs(xxt-xMarkerFlexE));
    %[~,indexFlexS] = min(abs(xxt+14));
    %[~,indexFlexE] = min(abs(xxt-14));

    clear marker_n xtV
    marker_n(:,:,:) = marker_cycle(indexFlexS+2:indexFlexE-2,:,:);
   
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
    clim1 = [-0.01 0.8 -5e-3 1.0 minForQ];
    clim2 = [ 0.01 1.2  5e-3 1.4 maxForQ];
    if(solid.thermoelasticFlag==1)
        nVar_xt_temp = nVar_xt;
    else
        nVar_xt_temp =  nVar_xt - 1;
    end
    figure;clf;
    tiledlayout(1,nVar_xt_temp,"TileSpacing","tight","Padding","compact");
    for iVar = 1:nVar_xt_temp
        nexttile;
        varValue = variableList_xt(iVar);
        clear xx xtV;
        xtV(:,:) = marker_n(:,varValue,:);
        if(varValue==4)
            nP_smooth = 3;
            KK = (1/nP_smooth^2)*ones(nP_smooth);
            
            for ii = 1:cycleGap
                xtV(:,ii) = smooth(xtV(:,ii),0.05,'loess');
            end

            for ii = 1:length(xtV(:,1))
                xtV(ii,:) = smooth(xtV(ii,:),0.05);
            end

            %xtV = conv2(xtV,KK,'same');
            %xtV = smooth(xtV,2); 
        end
        xx(:) = marker_n(:,1,1);
        [xxx,tt] = meshgrid(xx,t_cy);
        contourf(xxx,tt,xtV',200,"edgeColor","none")
        if(colormapType(iVar)==1)
            colormap(BlWtRd);
        else
            colormap("jet");
        end
        colorbar;
        hold on;
        for i=1:nLines
            yline(t_cy(snapshot_time_index(i)),"--k")
        end
     
        %ylim([0 300])
        caxis([clim1(iVar) clim2(iVar)])
        xlim([xMarkerFlexS xMarkerFlexE]/delta)
        if(iVar==1) 
            ylabel('$tU/a$','Interpreter','latex'); xlabel('$x/a$','Interpreter','latex')
        else
            yticks([])
        end
        set(gca, 'FontName', 'Times', 'ycolor','k',"FontSize",9)
        title(ttltxt{iVar},'Interpreter','latex')
        %saveas(gca,"user_xt_"+variableNames_xt(iVar)+"_"+folderNames{caseNumber},'fig')
    end
   
    set(gcf, 'Position',  [100, 500, 800, 250])
    saveas(gca,"cycle_avg_user_xt_"+variableNames_xt(iVar)+"_"+folderNames{caseNumber},'fig')

end
fclose(fidDampedSinosoid);
