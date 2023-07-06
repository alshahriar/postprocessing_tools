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
    panelLenght = 10;

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
        case_list_thermo_set_1;
    elseif(0)
        case_list_thermo_set_2;
    elseif(1)
        case_list_thermo_set_3;
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
nCases = length(folderNames);

defAtMiddleAllCases = cell(nCases,1);
def14thAllCases = cell(nCases,1);
def34thAllCases = cell(nCases,1);
freqValuesAllCases = cell(nCases,1);
peaksTimeIndexAllCases = cell(nCases,1);
valleysTimeIndexAllCases = cell(nCases,1);
coeffiecientsAllCases = zeros(nCases,5);
t_normAllCases = cell(nCases,1);
dampingRatioAllCases = cell(nCases,1);
nuAllCases = cell(nCases,1);
TAllCases = cell(nCases,1);

load(fullfile(dirName{1},folderNames{1},'marker.mat'));
multipleCasesPr = zeros(length(marker(:,1,1)),11,nCases);
clear tempData_1;

%%
averageDefatMiddle = zeros(1,nCases);
CauchyNumber = zeros(1,nCases);
freqValue = zeros(1,nCases);
%%
for caseNumber = 1:nCases
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
        plot(1:length(marker(midPointIndex,2,:)),squeeze(marker(midPointIndex,2,:)))
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

    [~,indexFlexS] = min(abs(xx-xMarkerFlexS));
    [~,indexFlexE] = min(abs(xx-xMarkerFlexE));


    if(solid.thermoelasticFlag==1)
        qxt = marker(indexFlexS:indexFlexE,11,:);
        nusslt = zeros(nTime,1);
        for i = 1:nTime
            nusslt(i) = mean(qxt(:,i));
        end
        nuAllCases{caseNumber} = nusslt;
    end

    Txt = marker(indexFlexS:indexFlexE,9,:);
    T = zeros(nTime,1);
    for i = 1:nTime
        T(i) = mean(Txt(:,i));
    end
    TAllCases{caseNumber} = T;

end
%%

if(set_number == 2)
    caseList = [2 3 4 5];
    nCasesToPlot = length(caseList);
    for i = 1:nCasesToPlot
        caseNuumber = caseList(i);
        disp([folderNames{caseNuumber}])
    end

    Legend{1} = "T-11a";
    Legend{2} = "T-12a";
    Legend{3} = "T-21a";
    Legend{4} = "T-22a";
    color = ['r';'g';'b';'k'];


    figure;
    nCasesToPlot = length(caseList);
    for i = 1:nCasesToPlot
        caseNuumber = caseList(i);
        nusslt  = nuAllCases{caseNuumber};
        t = t_normAllCases{caseNuumber};
        plot(t,smooth(nusslt,25),'Color',color(i)); hold on;
    end
    ylabel('$\langle Nu \rangle$','Interpreter','latex')
    xlabel('$tU/a$','Interpreter','latex')
    ylim([-0.0015 0.0015]*5e-3)
    grid on
    set(gca, 'FontName', 'Times','ycolor','k')
    set(gcf, 'Position',  [100, 500, 400, 250])
    xlim([0 300])
    legend(Legend,'Location','best')
    imageName = ['nussl_adiabatic_fixed'];
    box on; grid on;
    saveas(gcf,imageName,'fig')

    %% Temperature
    figure;
    nCasesToPlot = length(caseList);
    for i = 1:nCasesToPlot
        caseNuumber = caseList(i);
        T  = TAllCases{caseNuumber};
        t = t_normAllCases{caseNuumber};
        plot(t,T,'Color',color(i)); hold on
        grid on
        set(gca, 'FontName', 'Times','ycolor','k')
    end
    xlim([0 300])
    ylabel('$ \langle T \rangle /T_\infty$','Interpreter','latex')
    xlabel('$tU/a$','Interpreter','latex')
    %ylim([-0.003 0.003])
    grid on
    set(gca, 'FontName', 'Times','ycolor','k')
    set(gcf, 'Position',  [1300, 500, 400, 250])
    legend(Legend,'Location','best')
    imageName = ['T_adiabatic_fixed'];
    box on; grid on;
    saveas(gcf,imageName,'fig')

    %% Deflections - fixed - specific heat and expansion co
    figure;
    caseList = [2 3 4 5];
    nCasesToPlot = length(caseList);

    for i = 1:nCasesToPlot
        caseNuumber = caseList(i);
        def2  = defAtMiddleAllCases{caseNuumber};
        t = t_normAllCases{caseNuumber};
        plot(t,def2,'Color',color(i)); hold on;
    end

    for i = 1:nCasesToPlot
        caseNuumber = caseList(i);
        def2  = def34thAllCases{caseNuumber};
        t = t_normAllCases{caseNuumber};
        plot(t,def2,'Color',color(i)); hold on;
    end

    xlim([0 300])
    ylabel('$\Delta y /a$','Interpreter','latex')
    xlabel('$tU/a$','Interpreter','latex')
    %ylim([-0.0015 0.0015])
    grid on
    set(gca, 'FontName', 'Times','ycolor','k')
    set(gcf, 'Position',  [1300, 100, 400, 250])

    legend(Legend,'Location','best')
    imageName = ['diflection_fixed'];
    box on; grid on;
    saveas(gcf,imageName,'fig');

    %% Deflections - difference - fixed - specific heat and expansion co
    figure;
    minTimeIndex = 100000;
    for i = 1:nCasesToPlot
        caseNuumber = caseList(i);
        t = t_normAllCases{caseNuumber};
        minTimeIndex = min(minTimeIndex,length(t));
    end
    def = zeros(minTimeIndex,nCasesToPlot);
    ttrim = zeros(minTimeIndex,nCasesToPlot);
    for i = 1:nCasesToPlot
        caseNuumber = caseList(i);
        def(:,i) = def34thAllCases{caseNuumber}(1:minTimeIndex);
        ttrim(:,i) = t_normAllCases{caseNuumber}(1:minTimeIndex);
    end

    dd1 = def(:,1)-def(:,2);
    dd2 = def(:,1)-def(:,3);
    dd3 = def(:,1)-def(:,4);

    plot(ttrim(:,1),dd1,"r"); hold on;
    plot(ttrim(:,1),dd2,"b")
    plot(ttrim(:,1),dd3,"g")

    xlim([0 200])
    ylabel('$\Delta y /a$','Interpreter','latex')
    xlabel('$tU/a$','Interpreter','latex')
    %ylim([-0.0015 0.0015])
    grid on
    set(gca, 'FontName', 'Times','ycolor','k')
    set(gcf, 'Position',  [900, 100, 400, 250])
    legend("T11a-T12a","T11a-T21a","T11a-T22a",'Location','best')
    imageName = ['diflection_diff_fixed'];
    box on; grid on;
    saveas(gcf,imageName,'fig');


    %% Deflections - no-couple vs couple - fixed
    figure;
    caseList = [1 2];
    nCasesToPlot = length(caseList);

    for i = 1:nCasesToPlot
        caseNuumber = caseList(i);
        def2  = defAtMiddleAllCases{caseNuumber};
        t = t_normAllCases{caseNuumber};
        plot(t,def2,'Color',color(i)); hold on;
    end

    for i = 1:nCasesToPlot
        caseNuumber = caseList(i);
        def2  = def34thAllCases{caseNuumber};
        t = t_normAllCases{caseNuumber};
        plot(t,def2,'--','Color',color(i)); hold on;
    end

    xlim([0 300])
    ylabel('$y /a$','Interpreter','latex')
    xlabel('$tU/a$','Interpreter','latex')
    %ylim([-0.0015 0.0015])
    grid on
    set(gca, 'FontName', 'Times','ycolor','k')
    set(gcf, 'Position',  [1300, 100, 400, 250])
    clear Legend
    Legend{1} = "adiabatic"; Legend{2} = "coupled";
    legend(Legend,'Location','best')
    imageName = ['adiabaticVsCoupled_fixed'];
    box on; grid on;
    saveas(gcf,imageName,'fig');


    %% Deflections - no-couple vs couple - pinned
    figure;
    caseList = [6 7];
    nCasesToPlot = length(caseList);

    for i = 1:nCasesToPlot
        caseNuumber = caseList(i);
        def2  = defAtMiddleAllCases{caseNuumber};
        t = t_normAllCases{caseNuumber};
        plot(t,def2,'Color',color(i)); hold on;
    end

    for i = 1:nCasesToPlot
        caseNuumber = caseList(i);
        def2  = def34thAllCases{caseNuumber};
        t = t_normAllCases{caseNuumber};
        plot(t,def2,"--",'Color',color(i)); hold on;
    end

    xlim([0 300])
    ylabel('$y /a$','Interpreter','latex')
    xlabel('$tU/a$','Interpreter','latex')
    %ylim([-0.0015 0.0015])
    grid on
    set(gca, 'FontName', 'Times','ycolor','k')
    set(gcf, 'Position',  [900, 100, 400, 250])
    clear Legend
    Legend{1} = "adiabatic"; Legend{2} = "coupled";
    legend(Legend,'Location','best')
    imageName = ['adiabaticVsCoupled_pinned'];
    box on; grid on;
    saveas(gcf,imageName,'fig');

    %%
    caseList = 1:nCases;
    figure;
    nCasesToPlot = length(caseList);
    for i = 1:nCasesToPlot
        figure;
        caseNuumber = caseList(i);
        def  = defAtMiddleAllCases{caseNuumber};
        t = t_normAllCases{caseNuumber};
        plot(1:length(t),def);
        xlim([0 80])
        grid on
        set(gca, 'FontName', 'Times','ycolor','k')
        title("i"+caseNuumber)
    end
end

%%
%% Now set 3
if(set_number == 3)
    % Deflection
    figure;
    color = ['r';'g';'b';'k'];
    caseList = [1 2 3 4];
    clear Legend;   Legend{1} = "T-11a"; Legend{2} = "T-11c"; Legend{3} = "T-12c"; Legend{4} = "T-22c";
    nCasesToPlot = length(caseList);
    for i = 1:nCasesToPlot
        caseNuumber = caseList(i);
        def2  = defAtMiddleAllCases{caseNuumber};
        t = t_normAllCases{caseNuumber};
        plot(t,def2,'Color',color(i)); hold on;
    end
    xlim([0 300])
    ylabel('$y /a$','Interpreter','latex')
    xlabel('$tU/a$','Interpreter','latex')
    %ylim([-0.0015 0.0015])
    grid on
    set(gca, 'FontName', 'Times','ycolor','k')
    set(gcf, 'Position',  [1300, 100, 400, 250])
    legend(Legend,'Location','best')
    imageName = ['def_cold_vs_adb_fixed'];
    box on; grid on;
    saveas(gcf,imageName,'fig')

   figure;
    for i = 1:nCasesToPlot
        caseNuumber = caseList(i);
        def2  = def34thAllCases{caseNuumber};
        t = t_normAllCases{caseNuumber};
        plot(t,def2,"-",'Color',color(i)); hold on;
    end
    xlim([0 300])
    %ylim([-0.0015 0.0015])
    grid on
    ylabel('$y /a$','Interpreter','latex')
    xlabel('$tU/a$','Interpreter','latex')
    set(gca, 'FontName', 'Times','ycolor','k')
    set(gcf, 'Position',  [1300, 100, 400, 250])
    legend(Legend,'Location','best')
    imageName = ['def34_cold_vs_adb_fixed'];
    box on; grid on;
    saveas(gcf,imageName,'fig')

    %% Temp
    figure;
    color = ['r';'g';'b';'k'];
    caseList = [1 2 3 4];
    clear Legend;   Legend{1} = "T-11a"; Legend{2} = "T-11c"; Legend{3} = "T-12c"; Legend{4} = "T-22c";
    nCasesToPlot = length(caseList);
    for i = 1:nCasesToPlot
        caseNuumber = caseList(i);
        T2  = TAllCases{caseNuumber};
        t = t_normAllCases{caseNuumber};
        plot(t,T2,'Color',color(i)); hold on;
    end
    xlim([0 300])
    ylabel('$\langle T \rangle /T_\infty$','Interpreter','latex')
    xlabel('$tU/a$','Interpreter','latex')
    %ylim([-0.0015 0.0015])
    grid on
    set(gca, 'FontName', 'Times','ycolor','k')
    set(gcf, 'Position',  [900, 100, 400, 250])
    legend(Legend,'Location','best')

    axes("Position",[0.5 0.4 0.35 0.4]); box on;
    for i = 1:nCasesToPlot
        caseNuumber = caseList(i);
        T2  = TAllCases{caseNuumber};
        t = t_normAllCases{caseNuumber};
        plot(t,T2-1.05,'Color',color(i)); hold on;
    end
    xlim([150 300]); ylim([0.0027 0.00285])
    set(gca, 'FontName', 'Times','ycolor','k')

    imageName = ['temp_cold_vs_adb_fixed'];
    box on; grid on;
    saveas(gcf,imageName,'fig')

 %% Q
    figure;
    color = ['g';'b';'k'];
    caseList = [2 3 4];
    clear Legend;   Legend{1} = "T-11c"; Legend{2} = "T-12c"; Legend{3} = "T-22c";
    nCasesToPlot = length(caseList);
    for i = 1:nCasesToPlot
        caseNuumber = caseList(i);
        nu2  = nuAllCases{caseNuumber};
        t = t_normAllCases{caseNuumber};
        plot(t,nu2,'Color',color(i)); hold on;
    end
    xlim([0 300])
    ylabel('$\langle Nu \rangle$','Interpreter','latex')
    xlabel('$tU/a$','Interpreter','latex')
    %ylim([-0.0015 0.0015])
    grid on
    set(gca, 'FontName', 'Times','ycolor','k')
    set(gcf, 'Position',  [900, 100, 400, 250])
    legend(Legend,'Location','best')
    imageName = ['nu_cold_vs_adb_fixed'];
    box on; grid on;
    saveas(gcf,imageName,'fig')

    %% Heated vs adiabatic

    % def
    figure;
    color = ['r';'g'];
    caseList = [1 5];
    clear Legend;   Legend{1} = "T-11a"; Legend{2} = "T-11h";
    nCasesToPlot = length(caseList);
    for i = 1:nCasesToPlot
        caseNuumber = caseList(i);
        def2  = defAtMiddleAllCases{caseNuumber};
        t = t_normAllCases{caseNuumber};
        plot(t,def2,'Color',color(i)); hold on;
    end
    xlim([0 300])
    ylabel('$y /a$','Interpreter','latex')
    xlabel('$tU/a$','Interpreter','latex')
    %ylim([-0.0015 0.0015])
    grid on
    set(gca, 'FontName', 'Times','ycolor','k')
    set(gcf, 'Position',  [1300, 100, 400, 250])
    legend(Legend,'Location','best')
    imageName = ['def_hot_vs_adb_fixed'];
    box on; grid on;
    saveas(gcf,imageName,'fig')

    %% T Nu
    figure;
    nCasesToPlot = length(caseList);
    for i = 1:nCasesToPlot
        caseNuumber = caseList(i);
        T2  = TAllCases{caseNuumber};
        t = t_normAllCases{caseNuumber};
        plot(t,T2,"--",'Color',color(i)); hold on;
    end
    ylabel('$\langle T \rangle/ T_\infty$','Interpreter','latex')


    yyaxis right;
    for i = 1:nCasesToPlot
        caseNuumber = caseList(i);
        nu2  = nuAllCases{caseNuumber};
        t = t_normAllCases{caseNuumber};
        plot(t,nu2,"-",'Color',color(i)); hold on;
    end

    xlim([0 300])
    ylabel('$\langle Nu \rangle$','Interpreter','latex')
    xlabel('$tU/a$','Interpreter','latex')
    %ylim([-0.0015 0.0015])
    grid on
    set(gca, 'FontName', 'Times','ycolor','k')
    set(gcf, 'Position',  [900, 100, 400, 250])
    legend(Legend,'Location','best')
    imageName = ['nu_T_hot_vs_adb_fixed'];
    box on; grid on;
    saveas(gcf,imageName,'fig')

end