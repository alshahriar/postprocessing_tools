clear all
dirLocT = "/gpfs/home/as17r/thermo_movie/33_2/temp";
dirLocV = "/gpfs/home/as17r/thermo_movie/33_2/vort";
dirLocOut = "/gpfs/home/as17r/thermo_movie/33_2";

fname = "visit*.png";

U = 2;
fluid.ntec = 500;
fluid.dt = 0.002;
a = 10;

fullname = fullfile(dirLocT,fname);
filesInfoT = dir(fullname);
nImagesT = length(filesInfoT);

fullname = fullfile(dirLocV,fname);
filesInfoV = dir(fullname);
nImagesV = length(filesInfoV);
nImages = min([nImagesV,nImagesT]);
fig = figure; hold on;
caseNumber = 1;
dirName{caseNumber} = '/gpfs/home/as17r/engineering/ParametricStudy/AllCases2/Thermal_Coupling_Test_Cases/OldSetup/';
folderNames{caseNumber} = "test33_2_rerun";
load(fullfile(dirName{caseNumber},folderNames{caseNumber},'marker'),'marker');

nImages = 120;
t_norm = fluid.dt*fluid.ntec*U/a*(1:nImages);
xx = marker(:,1,1);
[~,midIndex] = min(abs(xx));
[~,three4Index] = min(abs(xx-2.5));
deflection(:,1) = (marker(midIndex,2,40:nImages+40)-1)/a;
deflection2(:,1) = (marker(three4Index,2,40:nImages+40)-1)/a;
yMin = min([deflection;deflection2]);
yMax = max([deflection;deflection2]);
stride = 5;

fullnameOutput = fullfile(dirLocOut,"combined");
outputVideo = VideoWriter(fullnameOutput,"Uncompressed AVI");
open(outputVideo)
for idx = 1:stride:nImages
    imT = imread(fullfile(filesInfoT(idx).folder,filesInfoT(idx).name));
    imV = imread(fullfile(filesInfoV(idx).folder,filesInfoV(idx).name));
    if(1)
        left = 415;
        right = 375;
        bottom = 180;
        top = 60;
        imT = imT(top:end-bottom,left:end-right,:);
        imV = imV(top:end-bottom,left:end-right,:);
        im = cat(2,imV,imT);
        figure(1); clf;
        plot(t_norm(1:idx),deflection(1:idx)); hold on;
        plot(t_norm(1:idx),deflection2(1:idx));
        xlim([t_norm(1) t_norm(nImages)])
        ylim([yMin yMax]);
        set(gcf,"Position",[100 100 length(im(1,:,1)) floor(length(im(:,1,1))/4)])
        ylabel('$y/a$','interpreter','latex')
        xlabel('$tU/a$','Interpreter','latex')
        legend("mid","3/4a","Location","northeast");
        F = getframe(gcf);
        [imd] = frame2im(F);
        im = cat(1,im,imd);
        im = imresize(im,0.5);
        %[A,map] = rgb2ind(im,256);
        writeVideo(outputVideo,im);
    end
end
close(outputVideo)
close;