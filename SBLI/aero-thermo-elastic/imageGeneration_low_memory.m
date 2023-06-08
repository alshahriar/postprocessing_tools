clear all
dirLocT = "/gpfs/home/as17r/thermo_movie/33_2/temp";
fname = "visit*.png";

fullname = fullfile(dirLoc,fname);
filesInfo = dir(fullname);
nImages = length(filesInfo);
fullnameOutput = fullfile(dirLoc,"movie.gif");
fig = figure; hold on;
for idx = 1:500
    im = imread(fullfile(filesInfo(idx).folder,filesInfo(idx).name));
    if(0)
        imshow(im,'InitialMagnification',100);
        if(idx==1)
            exportgraphics(gcf,fullnameOutput,'Append',false);
        else
            exportgraphics(gcf,fullnameOutput,'Append',true);
        end
    else
        left = 415;
        right = 375;
        bottom = 180;
        top = 60;
        im = im(top:end-bottom,left:end-right,:);
        im = imresize(im,0.5);
        [A,map] = rgb2ind(im,256);
        if idx == 1
            imwrite(A,map,fullnameOutput,"gif","LoopCount",Inf,"DelayTime",0.1);
        else
            imwrite(A,map,fullnameOutput,"gif","WriteMode","append","DelayTime",0.1);
        end
    end
end
close;