userDirCaseNumber = userDirCaseNumber + 1;
dirName{userDirCaseNumber} = '/gpfs/home/as17r/engineering/ParametricStudy/AllCases2/Thermal_Coupling_Test_Cases/OldSetup/';
folderNames{userDirCaseNumber} = "test29_corr";
iterS(userDirCaseNumber) = 1152; iterE(userDirCaseNumber) = 0;

userDirCaseNumber = userDirCaseNumber + 1;
dirName{userDirCaseNumber} = '/gpfs/home/as17r/engineering/ParametricStudy/AllCases2/Thermal_Coupling_Test_Cases/OldSetup/';
folderNames{userDirCaseNumber} = 'test32_corr'; %issue
iterS(userDirCaseNumber) = 448; iterE(userDirCaseNumber) = 0;

userDirCaseNumber = userDirCaseNumber + 1;
dirName{userDirCaseNumber} = '/gpfs/home/as17r/engineering/ParametricStudy/AllCases2/Thermal_Coupling_Test_Cases/OldSetup/';
folderNames{userDirCaseNumber} = 'test33_corr'; %issue
iterS(userDirCaseNumber) = 340; iterE(userDirCaseNumber) = 0;

userDirCaseNumber = userDirCaseNumber + 1;
dirName{userDirCaseNumber} = '/gpfs/home/as17r/engineering/ParametricStudy/AllCases2/Thermal_Coupling_Test_Cases/OldSetup/';
folderNames{userDirCaseNumber} = 'test32_corr2'; %fine
iterS(userDirCaseNumber) = 1; iterE(userDirCaseNumber) = 0;

userDirCaseNumber = userDirCaseNumber + 1;
dirName{userDirCaseNumber} = '/gpfs/home/as17r/engineering/ParametricStudy/AllCases2/Thermal_Coupling_Test_Cases/OldSetup/';
folderNames{userDirCaseNumber} = 'test33_corr2'; %fine
iterS(userDirCaseNumber) = 1; iterE(userDirCaseNumber) = 0;


userDirCaseNumber = userDirCaseNumber + 1;
dirName{userDirCaseNumber} = '/gpfs/home/as17r/engineering/ParametricStudy/AllCases2/Thermal_Coupling_Test_Cases/OldSetup/';
folderNames{userDirCaseNumber} = 'test29'; % fine
iterS(userDirCaseNumber) = 1; iterE(userDirCaseNumber) = 0;


userDirCaseNumber = userDirCaseNumber + 1;
dirName{userDirCaseNumber} = '/gpfs/home/as17r/engineering/ParametricStudy/AllCases2/Thermal_Coupling_Test_Cases/OldSetup/';
folderNames{userDirCaseNumber} = 'test32'; % fine
iterS(userDirCaseNumber) = 1; iterE(userDirCaseNumber) = 0;


userDirCaseNumber = userDirCaseNumber + 1;
dirName{userDirCaseNumber} = '/gpfs/home/as17r/engineering/ParametricStudy/AllCases2/Thermal_Coupling_Test_Cases/OldSetup/';
folderNames{userDirCaseNumber} = 'test33'; % fine
iterS(userDirCaseNumber) = 1; iterE(userDirCaseNumber) = 0;


userDirCaseNumber = userDirCaseNumber + 1;
dirName{userDirCaseNumber} = '/gpfs/home/as17r/engineering/ParametricStudy/AllCases2/Thermal_Coupling_Test_Cases/OldSetup/';
folderNames{userDirCaseNumber} = 'test32_2';
iterS(userDirCaseNumber) = 1; iterE(userDirCaseNumber) = 0;

userDirCaseNumber = userDirCaseNumber + 1;
dirName{userDirCaseNumber} = '/gpfs/home/as17r/engineering/ParametricStudy/AllCases2/Thermal_Coupling_Test_Cases/OldSetup/';
folderNames{userDirCaseNumber} = 'test33_2';
iterS(userDirCaseNumber) = 1; iterE(userDirCaseNumber) = 0;

userDirCaseNumber = userDirCaseNumber + 1;
dirName{userDirCaseNumber} = '/gpfs/home/as17r/engineering/ParametricStudy/AllCases2/Thermal_Coupling_Test_Cases/OldSetup/';
folderNames{userDirCaseNumber} = 'test29a';
iterS(userDirCaseNumber) = 1; iterE(userDirCaseNumber) = 0;

userDirCaseNumber = userDirCaseNumber + 1;
dirName{userDirCaseNumber} = '/gpfs/home/as17r/engineering/ParametricStudy/AllCases2/Thermal_Coupling_Test_Cases/OldSetup/';
folderNames{userDirCaseNumber} = 'test29b';
iterS(userDirCaseNumber) = 1; iterE(userDirCaseNumber) = 0;

for i =1:userDirCaseNumber
    uSubfolders{i} = 'flow_output';
    nvarList(i) = 11;
end
srow = 3056; erow = 4048;
nvarList(6) = 10;
nvarList(11) = 10;
nvarList(12) = 10;
nVar = 11;