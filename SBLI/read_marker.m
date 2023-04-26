clear all
%%
if(isunix)
    addpath('/gpfs/research/engineering/as17r/SBLI_Subroutines_Functions')
else
    addpath('SBLI_Subroutines_Functions')
end

ntec =  1000; % gap between each file
nVar = 10;

% smaller panel
srow = 3056;
erow = 4048;

%longer panel
srow = 2914;
erow = 3884;

% pannel fluttering case
srow = 4+1260;
erow = 4+2519;

% pannel fluttering - after static verification
srow = 4004;
erow = 8000;

nPoints = erow-srow+1;

userDefinedFolderOnly = 1;
userDirCaseNumber = 0;
if(userDefinedFolderOnly == 1)
    userDirCaseNumber = userDirCaseNumber + 1;
    uDirName{userDirCaseNumber} = '/gpfs/research/engineering/as17r/ParametricStudy/AllCases2/latestVersion/';
    uNameOfFolders{userDirCaseNumber} = 'newTahoeRigid28i';
    uSubfolders{userDirCaseNumber} = 'flow_output';
    %     userDirCaseNumber = userDirCaseNumber + 1;
    %     uDirName{userDirCaseNumber} = 'D:\CCNS\parametricStudy\AllCases2';
    %     uNameOfFolders{userDirCaseNumber} = 'H_T0600_E05000_veryLongRun_density_100';
    %
    %     userDirCaseNumber = userDirCaseNumber + 1;
    %     uDirName{userDirCaseNumber} = 'D:\CCNS\parametricStudy\AllCases2';
    %     uNameOfFolders{userDirCaseNumber} = 'H_T0600_E05000_veryLongRun_density_500_oldcode';
    nFolder = userDirCaseNumber;
else
    [dirName,folderNames,~,~] = getFileAndFolderName(1,0);
    nFolder = length(folderNames);
end

for caseNumber = 1: nFolder
    if(userDefinedFolderOnly)
        fullDirName{caseNumber} = fullfile(uDirName{caseNumber},uNameOfFolders{caseNumber},uSubfolders{caseNumber});
    else
        if(isunix)
            fullDirName{caseNumber} = fullfile('/gpfs/research/engineering/as17r/ParametricStudy', ...
                folderNames{caseNumber},uSubfolders{caseNumber});
        else
            fullDirName{caseNumber} = fullfile('E:\CCNS\ParametricStudy',folderNames{caseNumber},uSubfolders{caseNumber});
        end
    end
    fname = fullfile(fullDirName{caseNumber},'marker.');
    
    disp(['Data will be extracted from: ',fullDirName{caseNumber}])
    
    info = dir([fname,'*.dat']);
    if(isempty(info))
        disp('Err:--> Something is wrong - most probably the directory name');
        break;
    end
    
    disp(['Extracting Data from : ',info(1).folder])
    nMarkerFiles = length(info);
    disp(['Total number of files : ',int2str(nMarkerFiles)])
    n = 0;
    clear marker;
    r1 = srow; c1=1; r2=erow; c2=nVar;
    marker = zeros(nPoints,nVar,nMarkerFiles);
    for i=1:nMarkerFiles
        fullName = fullfile(info(i).folder,info(i).name);
        marker(:,:,i) = readmatrix(fullName, 'Range',[r1 c1 r2 c2]);
        disp(fullName);
    end
    save(fullfile(uDirName{caseNumber},uNameOfFolders{caseNumber},'marker'),'marker');
end