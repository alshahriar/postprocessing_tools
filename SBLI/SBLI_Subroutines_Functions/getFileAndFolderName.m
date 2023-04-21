function [dirName,folderNames,n_T_Case,n_E_Case] = getFileAndFolderName(n_T,n_E,dirNameParametricStudy,userDefinedDirectoryFlag,userDefineDirectoryAppendFlag,uDirName,uNameOfFolders,un_T,un_E)

nT = length(n_T);
nE = length(n_E);
caseNumber = 0;

if(userDefinedDirectoryFlag)
    if(userDefineDirectoryAppendFlag)
        for i = 1:nT
            for j = 1:nE
                caseNumber = caseNumber+1;
                folderNames{caseNumber} = sprintf('H_T%04d_E%05d',1000*n_T(i),n_E(j));
                n_T_Case(caseNumber) = n_T(i);
                n_E_Case(caseNumber) = n_E(j);
                dirName{caseNumber} = dirNameParametricStudy;

            end
        end
    end
    for i = 1: length(uNameOfFolders)
        caseNumber = caseNumber + 1;
        folderNames{caseNumber} = uNameOfFolders{i};
        n_T_Case(caseNumber) = un_T(i);
        n_E_Case(caseNumber) = un_E(i);
        dirName{caseNumber} = uDirName{i};
    end
else
    for i = 1:nT
        for j = 1:nE
            caseNumber = caseNumber+1;
            folderNames{caseNumber} = sprintf('H_T%04d_E%05d',1000*n_T(i),n_E(j));
            n_T_Case(caseNumber) = n_T(i);
            n_E_Case(caseNumber) = n_E(j);
            dirName{caseNumber} = dirNameParametricStudy;
        end
    end
end

end