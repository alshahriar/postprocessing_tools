% Reading the binary output of CCNS

clear all; close all; fclose all;

if(isunix)
    dirNameParametricStudy = pwd;
    folderName = '';
else
    dirNameParametricStudy = "D:\CCNS\";
    folderName = '';
end
%%

status = [0 0 0 0 1      0 0 0 0 0      1 1 0 0 0      0 1 0 0 1      1 1 1 0 0      0 0 0 0 0];

for caseNumber = 1:30
    if(status(caseNumber))
        folderName = sprintf("Case%02d",caseNumber);
        fname = fullfile(dirNameParametricStudy,folderName,'slice_xy_');
        disp(['Data will be extracted from: ',fname])
        info = dir(fname+"*.plt");
        if(isempty(info))
            disp('Err:--> Something is wrong - most probably the directory name');
            break;
        end
        % info = dir('fc*.plt');
        nFiles = length(info);
        fileName = info(1).name;
        fullFileName = fullfile(dirNameParametricStudy,folderName,fileName);
        
        [NI,NJ,NK,NUMVAR,VarName,~] = readHeaderInfo(fullFileName,1);
        
        xc = zeros(NI,NJ,1);
        yc = zeros(NI,NJ,1);
        zc = zeros(NI,NJ,1);
        rho = zeros(NI,NJ,nFiles);
        u = zeros(NI,NJ,nFiles);
        v = zeros(NI,NJ,nFiles);
        w = zeros(NI,NJ,nFiles);
        p = zeros(NI,NJ,nFiles);
        
        t = zeros(nFiles,1);
        
        for iFile = 1:nFiles
            fileName = info(iFile).name;
            disp(fileName)
            timestep = str2num(fileName(end-10:end-4));
            t(iFile) = timestep;
            fullFileName = fullfile(dirNameParametricStudy,folderName,fileName);
            [NI,NJ,NK,NUMVAR,VarName,fid] = readHeaderInfo(fullFileName,0);
            DATAMARKER = fread(fid, 1, 'real*4');
            ZONEMARKER = fread(fid, 1, 'real*4');
            NUMRP = fread(fid, 1, 'bit32');
            DATATYPE = fread(fid, NUMVAR, 'bit32');
            
            for j = 1:NJ
                xc(1:NI,j) = fread(fid, NI, 'double');
            end
            for j = 1:NJ
                yc(1:NI,j) = fread(fid, NI, 'double');
            end
            % Flow variables
            for j = 1:NJ
                rho(1:NI,j,iFile) = fread(fid, NI, 'double');
            end
            for j = 1:NJ
                u(1:NI,j,iFile) = fread(fid, NI, 'double');
            end
            for j = 1:NJ
                v(1:NI,j,iFile) = fread(fid, NI, 'double');
            end
            for j = 1:NJ
                w(1:NI,j,iFile) = fread(fid, NI, 'double');
            end
            for j = 1:NJ
                p(1:NI,j,iFile) = fread(fid, NI, 'double');
            end
            fclose(fid);
            if(0)
                figure(1); clf;
                contourf(xc(:,:,1),yc(:,:,1),p(:,:,1),50,'edgecolor','none');
                colormap('jet')
                caxis([0.4 2.1])
                colorbar
                axis equal
                xlabel('x/\delta')
                ylabel('y/\delta')
                set(gca, 'FontName', 'Times','FontSize',15);
                set(gcf, 'Position',  [100, 100, 1166, 300])
                saveas(gcf,sprintf('rho_%07d',timestep),'png')
            end
        end
        if(1)
            xx = xc(:,1);
            yy = yc(1,:);
            [~,js] = min(abs(yy+0.0));
            js = js - 0;
            [~,je] = min(abs(yy-10.0));
            is = 1;
            [~,ie] = min(abs(xx-30.0));
            xc = xc(is:ie,js:je);
            yc = yc(is:ie,js:je);
            rho = rho(is:ie,js:je,:);
            u = u(is:ie,js:je,:);
            v = v(is:ie,js:je,:);
            w = w(is:ie,js:je,:);
            p = p(is:ie,js:je,:);            
        end
        save("flowField_"+caseNumber+".mat",'xc','yc','rho','u','v','w','p','nFiles','t','-v7.3')
        fid = fopen('flowdata.bin',w); fwrite(fid,rho,'double'); fclose(fid);
    end
end

%%
function [NI,NJ,NK,NUMVAR,VarName,fid] = readHeaderInfo(fullFileName,closeFlag)
fid=fopen(fullFileName, 'rb');        % Open the file.
magicNumber = fread(fid, 8, 'char');  % Reading the magic number
%fprintf('%s',char(magicNumber))
indetifier = fread(fid, 1, 'bit32');  % Reading indentifier
% Title   = 'CCNS Car3D Combined#'
temp1 = fread(fid, 19, 'bit32');  % Reading title
%fprintf('%s',char(temp1))
Title = sprintf('%s',char(temp1));
%TEMP
temp2 = fread(fid, 1, 'bit32'); % zero value
%NUMVAR
NUMVAR = fread(fid, 1, 'bit32');
% VarName = 'x,y,rho,u,v,p,pp,wz,GC,Dil,xmuR#'
% x y z rho u v w p T BL Dil xmuR
temp4 = fread(fid, 23, 'bit32');
%fprintf('%s',char(temp4))
VarName = sprintf('%s',char(temp4));
temp5 = fread(fid, 1, 'bit32'); % zero value
ZONEMARKER = fread(fid, 1, 'real*4');
temp = fread(fid, 6, 'bit32');  % Reading title
%fprintf('%s',char(temp))
ZONENAME = sprintf('%s',char(temp));
temp = fread(fid, 1, 'bit32'); % zero value
IFORMAT = fread(fid, 1, 'bit32');
ZONECOLOR = fread(fid, 1, 'bit32');
NI = fread(fid, 1, 'bit32');
NJ = fread(fid, 1, 'bit32');
NK = fread(fid, 1, 'bit32');
if(closeFlag==1)
    fclose(fid);
end
end

