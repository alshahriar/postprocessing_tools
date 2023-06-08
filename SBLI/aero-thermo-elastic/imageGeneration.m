clear all
dirLoc = "/gpfs/home/as17r/thermo_movie/33_2/temp";
fname = "movie*.png";

fullname = fullfile(dirLoc,fname);
filesInfo = dir(fullname);
nImages = length(filesInfo);

fig = figure;
for idx = 1:nImages
    im{idx} = imread(fullfile(filesInfo(idx).folder,filesInfo(idx).name));
end
close;

filename = fullfile(dirLoc,"movie.gif"); % Specify the output file name
for idx = 1:nImages
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(A,map,filename,"gif","LoopCount",Inf,"DelayTime",1);
    else
        imwrite(A,map,filename,"gif","WriteMode","append","DelayTime",1);
    end
end