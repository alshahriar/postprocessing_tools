set_number = 2;

% Fixed-Fixed
userDirCaseNumber = userDirCaseNumber + 1;
dirName{userDirCaseNumber} = '/gpfs/home/as17r/engineering/ParametricStudy/AllCases2/Thermal_Coupling_Test_Cases/OldSetup/';
folderNames{userDirCaseNumber} = 'test29_rerun'; % fine
iterS(userDirCaseNumber) = 40; iterE(userDirCaseNumber) = 0;

userDirCaseNumber = userDirCaseNumber + 1;
dirName{userDirCaseNumber} = '/gpfs/home/as17r/engineering/ParametricStudy/AllCases2/Thermal_Coupling_Test_Cases/OldSetup/';
folderNames{userDirCaseNumber} = 'test32_2_rerun';
iterS(userDirCaseNumber) = 40; iterE(userDirCaseNumber) = 0;

userDirCaseNumber = userDirCaseNumber + 1;
dirName{userDirCaseNumber} = '/gpfs/home/as17r/engineering/ParametricStudy/AllCases2/Thermal_Coupling_Test_Cases/OldSetup/';
folderNames{userDirCaseNumber} = 'test32_3_rerun';
iterS(userDirCaseNumber) = 20; iterE(userDirCaseNumber) = 0;

userDirCaseNumber = userDirCaseNumber + 1;
dirName{userDirCaseNumber} = '/gpfs/home/as17r/engineering/ParametricStudy/AllCases2/Thermal_Coupling_Test_Cases/OldSetup/';
folderNames{userDirCaseNumber} = 'test32_4_rerun';
iterS(userDirCaseNumber) = 1; iterE(userDirCaseNumber) = 0;

userDirCaseNumber = userDirCaseNumber + 1;
dirName{userDirCaseNumber} = '/gpfs/home/as17r/engineering/ParametricStudy/AllCases2/Thermal_Coupling_Test_Cases/OldSetup/';
folderNames{userDirCaseNumber} = 'test32_rerun'; % fine
iterS(userDirCaseNumber) = 40; iterE(userDirCaseNumber) = 0;

% Pinned-Pinned
userDirCaseNumber = userDirCaseNumber + 1;
dirName{userDirCaseNumber} = '/gpfs/home/as17r/engineering/ParametricStudy/AllCases2/Thermal_Coupling_Test_Cases/OldSetup/';
folderNames{userDirCaseNumber} = "test29_corr_rerun";
iterS(userDirCaseNumber) = 40; iterE(userDirCaseNumber) = 0;

userDirCaseNumber = userDirCaseNumber + 1;
dirName{userDirCaseNumber} = '/gpfs/home/as17r/engineering/ParametricStudy/AllCases2/Thermal_Coupling_Test_Cases/OldSetup/';
folderNames{userDirCaseNumber} = 'test32_corr2_rerun'; %issue
iterS(userDirCaseNumber) = 40; iterE(userDirCaseNumber) = 0;

userDirCaseNumber = userDirCaseNumber + 1;
dirName{userDirCaseNumber} = '/gpfs/home/as17r/engineering/ParametricStudy/AllCases2/Thermal_Coupling_Test_Cases/OldSetup/';
folderNames{userDirCaseNumber} = 'test32_corr3_rerun'; %fine
iterS(userDirCaseNumber) = 1; iterE(userDirCaseNumber) = 0;

userDirCaseNumber = userDirCaseNumber + 1;
dirName{userDirCaseNumber} = '/gpfs/home/as17r/engineering/ParametricStudy/AllCases2/Thermal_Coupling_Test_Cases/OldSetup/';
folderNames{userDirCaseNumber} = 'test32_corr4_rerun'; %issue
iterS(userDirCaseNumber) = 1; iterE(userDirCaseNumber) = 0;

for i =1:userDirCaseNumber
    uSubfolders{i} = 'flow_output';
    nvarList(i) = 11;
end
srow = 3056; erow = 4048;