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
        case_list_thermo_set_1;
    elseif(0)
        case_list_thermo_set_2;
    elseif(1)
        case_list_thermo_set_3;
    end
else
    case_list_old_study_set_1;
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
        plot(1:length(marker(midPointIndex,2,:)),squeeze(marker(midPointIndex,2,:)))
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
  
    %% xt plot of variables
    variableList_xt = [9 11];
    variableNames_xt = 'yvPTQ';

    colormapType = [1 1 2 2 2];

    J = customcolormap([0 0.5 1], {'#0000FF','#FFFFFF','#FF0000'},100);
    BlWtGr = colormap(J);

    nVar_xt = length(variableList_xt);
    xxt(:) = marker(:,1,1);
%      [~,indexFlexS] = min(abs(xxt-xMarkerFlexS));
%      [~,indexFlexE] = min(abs(xxt-xMarkerFlexE));

     [~,indexFlexS] = min(abs(xxt+14));
     [~,indexFlexE] = min(abs(xxt-14));

    clear marker_n xtV
    marker_n(:,:,:) = marker(indexFlexS+2:indexFlexE-2,:,:);
    marker_n(:,2,:) =  marker_n(:,2,:)-1;

    if(solid.thermoelasticFlag==1)
        xtV(:,:) = marker_n(:,11,:);
        %minForQ = min(min(xtV));
        %maxForQ = max(max(xtV));
        minForQ = -5e-5;
        maxForQ = 5e-5;
    else
        minForQ = 1;
        maxForQ = 2;
    end

    xtV(:,:) = marker_n(:,9,:)-0;
    minForT = min(min(xtV));
    maxForT = max(max(xtV));
    clim1 = [-0.25 -5e-3 0.65 minForT minForQ];
    clim2 = [0.25 5e-3 1.1 maxForT maxForQ];
    if(solid.thermoelasticFlag==1)
        nVar_xt_temp = nVar_xt;
    else
        nVar_xt_temp =  nVar_xt - 1;
    end

    for iVar = 1:nVar_xt_temp
        figure(9982134);clf;
        varValue = variableList_xt(iVar);
        clear xx xtV;
        xtV(:,:) = marker_n(:,varValue,:);
        xx(:) = marker_n(:,1,1);
        [xxx,tt] = meshgrid(xx,t_norm);
        contourf(xxx,tt,xtV',500,"edgeColor","none")
        if(colormapType(iVar)==1)
            colormap(BlWtGr);
        else
            colormap("jet");
        end
        colorbar;
        set(gca, 'FontName', 'Times', 'ycolor','k')
        set(gcf, 'Position',  [100, 100, 400, 250])
        ylim([0 300])
        caxis([clim1(iVar) clim2(iVar)])
        %xlim([xMarkerFlexS xMarkerFlexE])
        saveas(gca,"user_xt_"+variableNames_xt(iVar)+"_"+folderNames{caseNumber},imageFormat)
        %saveas(gca,"user_xt_"+variableNames_xt(iVar)+"_"+folderNames{caseNumber},'fig')
    end

end
