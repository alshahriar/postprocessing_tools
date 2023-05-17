%% Data-driven Resolvent Analysis for Separation Bubble
% Tso-Kang Wang, 04/26/2021
clear
close all; clc

% Visualization parameters
fontname = 'CMU Serif';
set(0,'DefaultAxesFontName',fontname,'DefaultTextFontName',fontname);
set(0,'DefaultAxesFontSize', 14);

readinput = 1;

%% Input
% The snapshots should be structured as {x1, x2,..., xj}
% each xj is the reshape of the flowfield data of jth snapshot

% Input parameters
step = 0.02;
dframe = 50;
dt = step*dframe;

nx=504; ny=128;

flist = 20300:50:25000;
filecount = length(flist);

caselist = [1:5,7,9:12];
casecount = length(caselist);

%% Read flow data
% read grid
inputpath = 'G:\Compliant_Wall\DataDrivenResolvent\rigid_1\qSECtxt0020025.plt';
[x, y, ~, ~, ~, ~] = flowData(inputpath, nx, ny);
xgrid=reshape(x,[nx, ny]); ygrid=reshape(y,[nx, ny]);

% read all input
if readinput==1    
    for icase = 1:casecount
        disp(['flow data_',num2str(icase),'/',num2str(casecount)])
        
        for ifile = 1:filecount            
            idata = (icase-1)*filecount + ifile;
            filepath = strcat('G:\Compliant_Wall\DataDrivenResolvent\rigid_', ...
                sprintf('%d',(caselist(icase))),'\',sprintf('qSECtxt%07d.plt', flist(ifile)));
            [~, ~, U(:,idata), V(:,idata), Vor(:,idata)] = flowData(filepath, nx,ny);           
        end

    end        
    disp('read data DONE')    
end

% create two sets of data matrices {x1, x2,..., xj-1} and {x2, x3,...,xj}
for icase = 1:casecount
    vor1(:,(icase-1)*(filecount-1)+1: icase*(filecount-1)) = Vor(:,(icase-1)*filecount+1: icase*filecount-1);
    vor2(:,(icase-1)*(filecount-1)+1: icase*(filecount-1)) = Vor(:,(icase-1)*filecount+2: icase*filecount);
end

%% DMD modes

% Solve the SVD problem
[U,Sigma,V] = svd(vor1,0);

% Calculate companion matrix
S = U'*vor2*V*diag(1./diag(Sigma));

% Compute DMD modes with companion matrix
[eigenvec,lambda] = eig(S);
phi = U*eigenvec;

% magnitude of the DMD modes
magnitude = pinv(phi)*vor1(:,1);
mag = abs(magnitude);
[mag_sort, sort_index] = sort(mag, 'descend'); % get the index of the largest modes

% frequencies and growth/decay rates
freq = angle(diag(lambda))/(2*pi*dt);
growth = log(abs(diag(lambda)))/dt;

for k = 1:4
    
    figure(k)
    maxplot = mean(abs(phi(:,sort_index(k)))) + std(abs(phi(:,sort_index(k))));
    contourf(xgrid,ygrid,reshape(real(phi(:,sort_index(k))),[nx,ny]),512,'LineStyle','none')
%     colormap(b2r(-maxplot,maxplot)), colorbar
    axis equal
    
    title(strcat('DMD mode ',num2str(k)))
     
end

%% DMD
[lambda, V, W, b] = DMD(vor1, vor2, dt, red);

% figure(1)
% contourf(x_Cart,y_Cart,reshape(real(V(:,1)),[ns,nr_max]), 512,'LineStyle','none')
% axis equal
% colorbar; caxis([-0.01, 0.01]);
% xlim([-10,2]);
% ylim([-1,3]);

%% Resolvent Gains
% Q = nthroot(Jcb, 2)'*nthroot(Jcb, 2);
% [F, flag] = chol(Q);
% 
% F = eye(size(vor2,2));

Q = V'*V;
[F, flag] = chol(Q);
Fi = inv(F);

wspan = 0.02:0.02:3;
lgain = 4;
Lmb = diag(lambda);

R = zeros([lgain, length(wspan)]);
for i = 1:length(wspan)
   
    i
    Rtmp = svd(F*inv(-1i*wspan(i)*eye(length(lambda))-Lmb)*Fi);
    R(:,i) = Rtmp(1:lgain).^2;
    
end

figure(5)
plot(wspan, R)
set(gca, 'YScale', 'log')
% ylim([5e2, 5e4])
title('Resolvent Gain', 'interpreter','latex')
xlabel('$\omega$','interpreter', 'latex')
ylabel('$\sigma_j$','interpreter', 'latex')
saveas(gcf, 'resolvent_gain.png')

%% Optimal forcing modes
omg = 1.28;
fred = 24;

% [PsiF, sigmaF, PhiF] = opt_forcing(lambda(1:fred), real(phi(:,1:fred)), eye(1), omg);
[PsiF, sigmaF, PhiF] = opt_forcing(lambda(1:fred), V(:,1:fred), eye(1), omg);

% [PsiF, sigmaF, PhiF] = svd(F*inv(-1i*omg*eye(length(lambda))-Lmb)*Fi, 0);
% 
% PsiF = V*Fi*PsiF;
% PhiF = V*Fi*PhiF;

% forcing mode
for k = 1:4
    figure(6+(k-1))
    contourf(x_Cart,y_Cart,reshape(real(PhiF(:,k)),[ns,nr_max]), 64,'LineStyle','none')
    colorbar; colormap(b2r(-0.02, 0.02));
    axis equal
    xlim([-10,2]);
    ylim([-1,3]);
    title(sprintf('forcing mode %d, f %02f', k, omg))
    saveas(gcf, sprintf('forcing_%d_f_%02f.png', k, omg))
end

% response mode
for k =1:4
    figure(10+(k-1))
    contourf(x_Cart,y_Cart,reshape(real(PsiF(:,k)),[ns,nr_max]), 64,'LineStyle','none')
    colorbar; colormap(b2r(-0.02, 0.02));
    axis equal
    xlim([-10,2]);
    ylim([-1,3]);
    title(sprintf('response mode %d, f %f', k, omg))
    saveas(gcf, sprintf('response_%d_f_%02f.png', k, omg))
end
