%readWallPressure
clear all; close all; fclose all;

hpcFlag = 0;

if(hpcFlag==1)
    dirNameParametricStudy = pwd;
else
%     dirNameParametricStudy = 'C:\Users\alsha\OneDrive - Florida State University\swbli\';
%     folderName = 'readingCCNSBinaryOutput';

    dirNameParametricStudy = 'G:\DataDrivenResolvent\3D_SBLI\randominput\';
    folderName = '1\runSTBLI';
end
%%

for caseNumber = 1: 1

    fname = fullfile(dirNameParametricStudy,folderName,'f_flow');

    disp(['Data will be extracted from: ',fname])
    
    info = dir([fname,'*.plt']);
    
    if(isempty(info))
        disp('        Err:--> Something is wrong - most probably the directory name');
        break;
    end
    
    % info = dir('fc*.plt');
    nFiles = length(info);
    
    fileName = info(1).name;
    fullFileName = fullfile(dirNameParametricStudy,folderName,fileName);
    fid=fopen(fullFileName, 'rb');        % Open the file.
    magicNumber = fread(fid, 8, 'char');  % Reading the magic number
    fprintf('%s',char(magicNumber))
    indetifier = fread(fid, 1, 'bit32');  % Reading indentifier
    % Title   = 'CCNS Car3D Combined#'
    temp1 = fread(fid, 19, 'bit32');  % Reading title
    fprintf('%s',char(temp1))
    Title = sprintf('%s',char(temp1));
    %TEMP
    temp2 = fread(fid, 1, 'bit32'); % zero value
    %NUMVAR
        NUMVAR = fread(fid, 1, 'bit32');
    % VarName = 'x,y,rho,u,v,p,pp,wz,GC,Dil,xmuR#'
    temp4 = fread(fid, 31, 'bit32'); 
    fprintf('%s',char(temp4))
    VarName = sprintf('%s',char(temp4));
    temp5 = fread(fid, 1, 'bit32'); % zero value
    ZONEMARKER = fread(fid, 1, 'real*4');
    temp = fread(fid, 6, 'bit32');  % Reading title
    fprintf('%s',char(temp))
    ZONENAME = sprintf('%s',char(temp));
    temp = fread(fid, 1, 'bit32'); % zero value    
    IFORMAT = fread(fid, 1, 'bit32');
    ZONECOLOR = fread(fid, 1, 'bit32');
    NI = fread(fid, 1, 'bit32');
    NJ = fread(fid, 1, 'bit32');
    NK = fread(fid, 1, 'bit32');
    fclose(fid);
    
    xc = zeros(NI,NJ,NK);
    yc = zeros(NI,NJ,NK);
    zc = zeros(NI,NJ,NK);    
    rho = zeros(NI,NJ,NK);
    u = zeros(NI,NJ,NK);
    v = zeros(NI,NJ,NK);
    w = zeros(NI,NJ,NK);
    p = zeros(NI,NJ,NK);
    T = zeros(NI,NJ,NK);
    GC = zeros(NI,NJ,NK);
    Dil = zeros(NI,NJ,NK);
    xmuR = zeros(NI,NJ,NK);
    
    rhoAll = zeros(NI,NJ,NK,nFiles);
    uAll = zeros(NI,NJ,NK,nFiles);
    vAll = zeros(NI,NJ,NK,nFiles);
    pAll = zeros(NI,NJ,NK,nFiles);
    gmapAll = zeros(NI,NJ,NK,nFiles);
    t = zeros(nFiles,1);
    
    for iFile = 1:nFiles
        fileName = info(iFile).name;
        disp(fileName)
        timestep = str2num(fileName(end-10:end-4));
        t(iFile) = timestep;
        fullFileName = fullfile(dirNameParametricStudy,folderName,fileName);
        fid=fopen(fullFileName, 'rb');        % Open the file.
        magicNumber = fread(fid, 8, 'char');  % Reading the magic number
        fprintf('%s',char(magicNumber))
        indetifier = fread(fid, 1, 'bit32');  % Reading indentifier
        % Title   = 'CCNS Car3D Combined#'
        temp = fread(fid, 19, 'bit32');  % Reading title
        fprintf('%s',char(temp))
        Title = sprintf('%s',char(temp));
        %TEMP
        temp = fread(fid, 1, 'bit32'); % zero value
        %NUMVAR
        NUMVAR = fread(fid, 1, 'bit32');
        % VarName = 'x,y,rho,u,v,p,pp,wz,GC,Dil,xmuR#'
        temp = fread(fid, 31, 'bit32');  % Reading title
        fprintf('%s',char(temp))
        VarName = sprintf('%s',char(temp));
        temp = fread(fid, 1, 'bit32'); % zero value
        ZONEMARKER = fread(fid, 1, 'real*4');
        temp = fread(fid, 6, 'bit32');  % Reading title
        fprintf('%s',char(temp))
        ZONENAME = sprintf('%s',char(temp));
        temp = fread(fid, 1, 'bit32'); % zero value
        IFORMAT = fread(fid, 1, 'bit32');
        ZONECOLOR = fread(fid, 1, 'bit32');
        NI = fread(fid, 1, 'bit32');
        NJ = fread(fid, 1, 'bit32');
        NK = fread(fid, 1, 'bit32');
        DATAMARKER = fread(fid, 1, 'real*4');
        ZONEMARKER = fread(fid, 1, 'real*4');
        NUMRP = fread(fid, 1, 'bit32');
        DATATYPE = fread(fid, NUMVAR, 'bit32');
        
        for k = 1:NK
            for j = 1:NJ
                xc(1:NI,j,k) = fread(fid, NI, 'double');
            end
        end
        for k = 1:NK
            for j = 1:NJ
                yc(1:NI,j,k) = fread(fid, NI, 'double');
            end
        end
        for k = 1:NK
            for j = 1:NJ
                zc(1:NI,j,k) = fread(fid, NI, 'double');
            end
        end
        for k = 1:NK
            for j = 1:NJ
                rho(1:NI,j,k) = fread(fid, NI, 'double');
            end
        end
        for k = 1:NK
            for j = 1:NJ
                u(1:NI,j,k) = fread(fid, NI, 'double');
            end
        end
        for k = 1:NK
            for j = 1:NJ
                v(1:NI,j,k) = fread(fid, NI, 'double');
            end
        end
        for k = 1:NK
            for j = 1:NJ
                w(1:NI,j,k) = fread(fid, NI, 'double');
            end
        end
        for k = 1:NK
            for j = 1:NJ
                p(1:NI,j,k) = fread(fid, NI, 'double');
            end
        end
        for k = 1:NK
            for j = 1:NJ
                T(1:NI,j,k) = fread(fid, NI, 'double');
            end
        end
        for k = 1:NK
            for j = 1:NJ
                GC(1:NI,j,k) = fread(fid, NI, 'double');
            end
        end
        for k = 1:NK
            for j = 1:NJ
                Dil(1:NI,j,k) = fread(fid, NI, 'double');
            end
        end
        for k = 1:NK
            for j = 1:NJ
                xmuR(1:NI,j,k) = fread(fid, NI, 'double');
            end
        end
        fclose(fid);
        
        rhoAll(:,:,:,iFile) = rho;
        uAll(:,:,:,iFile) = u;
        vAll(:,:,:,iFile) = v;
        pAll(:,:,:,iFile) = p;        
        gmapAll(:,:,:,iFile) = GC;
        
        
        figure(1); clf;
        contourf(xc(:,:,1),yc(:,:,1),rho(:,:,1),50,'edgecolor','none');
        colormap('jet')
        % caxis([0.4 2.1])
        colorbar
        axis equal
        xlabel('x/\delta')
        ylabel('y/\delta')
        set(gca, 'FontName', 'Times','FontSize',15);
        set(gcf, 'Position',  [100, 100, 1166, 300])
        saveas(gcf,sprintf('GMAP_%07d',timestep),'png')
    end    
    %save(['flowFieldData_',folderName],'xc','yc','rhoAll','uAll','vAll','pAll','t')
end





