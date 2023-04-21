function [tahoeParameters] = loadTahoeInput(dirName,foldername)
inPutFileName = 'input_tahoe.dat';
fullFileName = fullfile(dirName,foldername,inPutFileName);
data = getTheDataFromTheLine(fullFileName,3);
tahoeParameters.nstart_tahoe = data(1);
tahoeParameters.ndump = data(2);
tahoeParameters.patm = data(4);
tahoeParameters.pCavity = data(5);
data = getTheDataFromTheLine(fullFileName,25);
tahoeParameters.rho = data(1);
data = getTheDataFromTheLine(fullFileName,27);
tahoeParameters.E = data(1);
data = getTheDataFromTheLine(fullFileName,29);
tahoeParameters.nu = data(1);
end % loadCaseData function

function data = getTheDataFromTheLine(fullFileName,n)
fid = fopen(fullFileName);
for i=1:n
    tline = fgetl(fid);
end
data = extractNumFromStr(tline);
fclose(fid);
end

function numArray = extractNumFromStr(str)
str1 = regexprep(str,'[,;=]', ' ');
str2 = regexprep(regexprep(str1,'[^- 0-9.eE(,)/]',''), ' \D* ',' ');
str3 = regexprep(str2, {'\.\s','\E\s','\e\s','\s\E','\s\e'},' ');
numArray = str2num(str3);
end