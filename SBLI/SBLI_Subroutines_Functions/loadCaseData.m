function [simParameters] = loadCaseData(dirName,foldername)
inPutFileName = 'CCNS.inp';
fullFileName = fullfile(dirName,foldername,inPutFileName);
data = getTheDataFromTheLine(fullFileName,2);
simParameters.ntec = data(2);
simParameters.dt = data(3);
data = getTheDataFromTheLine(fullFileName,8);
simParameters.Uinf = data(1);
simParameters.Re = data(2);
data = getTheDataFromTheLine(fullFileName,32);
simParameters.delta = data(4);
simParameters.Twall = data(8);
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