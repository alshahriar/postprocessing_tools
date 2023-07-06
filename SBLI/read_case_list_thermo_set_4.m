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