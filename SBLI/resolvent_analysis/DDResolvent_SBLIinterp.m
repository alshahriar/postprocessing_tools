%% Data-driven Resolvent Analysis for 3D SBLI
% Tso-Kang Wang, 08/01/2022
clear
close all; clc

% Visualization parameters
font_name='CMU Serif';
font_axes=8;
font_legend=8; % if there are many or 8 if there are few
font_axestitle=8;
font_title=8;

marginleft=0.1;
marginright=0.1;
marginbottom=0.1;
margintop=0.1;

widthIn = 3.5;  %half page of JFM
heightIn = 2.6;    %1.8-2.6 is a nice range

% Switch
readinput = 1;
plotDMD = 1;

%% Input
% The snapshots should be structured as {x1, x2,..., xj}
% each xj is the reshape of the flowfield data of jth snapshot

% Input parameters
step = 0.02; % dt*u/delta
dframe = 50;
dt = step*dframe;

xs=-30; xe=30;  xres=510;   nx=180; % interpolation resolution
ys=0.01; ye=9;   yres=225;   ny=36;

flist = 17000:50:20300;
filecount = length(flist);

caselist = [1:18,20,21];
casecount = length(caselist);

red = 400; % reduction of dimension

% Declare variable size
rhoALL = zeros([nx*ny, filecount, casecount]);

%% Read flow data

% create grid
xx = linspace(xs, xe, nx);
yy = linspace(ys, ye, ny);
[xgrid, ygrid] = meshgrid(xx, yy);

% read all input
if readinput==1
    
    for icase = 1:casecount        
        disp(['flow data_',num2str(icase),'/',num2str(casecount)])
        
        load(sprintf('flowFieldData_%d.mat', caselist(icase)));        
        for ifile = 1:filecount
            xorg=squeeze(x2D(:,:,1,ifile)); yorg=squeeze(y2D(:,:,1,ifile));
            
            rhoorg=squeeze(rho2D(:,:,1,ifile));
            rho = interp2(xorg', yorg', rhoorg', xgrid, ygrid);
            rhoALL(:,ifile,icase) = reshape(rho', [nx*ny, 1]);            
        end        
        
    end
         
    disp('read data DONE')
end

% create two sets of data matrices {x1, x2,..., xj-1} and {x2, x3,...,xj}
for icase = 1:casecount
    rho1(:,(icase-1)*(filecount-1)+1: icase*(filecount-1)) = ...
        rhoALL(:,1:filecount-1,icase);
    rho2(:,(icase-1)*(filecount-1)+1: icase*(filecount-1)) = ...
        rhoALL(:,2:filecount,icase);
end


%% DMD modes
if plotDMD==1    
    % Solve the SVD problem
    [U,Sigma,V] = svd(rho1,0);
    
    % Calculate companion matrix
    S = U'*rho2*V*diag(1./diag(Sigma));
    
    % Compute DMD modes with companion matrix
    [eigenvec,lambda] = eig(S);
    phi = U*eigenvec;
    
    % magnitude of the DMD modes
    magnitude = pinv(phi)*rho1(:,1);
    mag = abs(magnitude);
    [mag_sort, sort_index] = sort(mag, 'descend'); % get the index of the largest modes
    
    % frequencies and growth/decay rates
    freq = angle(diag(lambda))/(2*pi*dt);
    growth = log(abs(diag(lambda)))/dt;
    
    for k = 1:4
        figure(k)
        maxplot = mean(abs(phi(:,sort_index(k)))) + std(abs(phi(:,sort_index(k))));
        contourf(xgrid',ygrid',reshape(real(phi(:,sort_index(k))),[nx,ny]),512,'LineStyle','none')
        %     colormap(b2r(-maxplot,maxplot)), colorbar
        axis equal
        
        title(strcat('DMD mode ',num2str(k)))
    end    
end

%% DMD
[lambda, V, W, b] = DMD(rho1, rho2, dt, red);

% figure(1)
% contourf(xgrid,ygrid,reshape(real(W(:,1)),[nx,ny]), 512,'LineStyle','none')
% axis equal
% colorbar; 

%% Resolvent Gains
% Q = nthroot(Jcb, 2)'*nthroot(Jcb, 2);
% [F, flag] = chol(Q);
% 
% F = eye(size(rho2,2));

Q = V'*V;
[F, flag] = chol(Q');
Fi = F\eye(size(F,1));

wspan = -1.2:0.01:1.2;
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
% title('Resolvent Gain', 'interpreter','latex')
set(gca,'Fontname',font_name,'Fontsize',font_axes,'Box','off','LineWidth',1)
xlim([-1.2, 1.2])
xlabel('$\omega$','Interpreter','Latex','fontsize',font_axes)
ylabel('$\sigma_j$','Interpreter','Latex','fontsize',font_axes)
plot_setting

print(gcf,'./figures/resolvent_gain.png','-dpng')

%% Optimal forcing modes
omg = 0.09;
fred = 400;

% [PsiF, sigmaF, PhiF] = opt_forcing(lambda(1:fred), real(phi(:,1:fred)), eye(1), omg);
[PsiF, sigmaF, PhiF] = opt_forcing(lambda(1:fred), V(:,1:fred), omg);

% [PsiF, sigmaF, PhiF] = svd(F*inv(-1i*omg*eye(length(lambda))-Lmb)*Fi, 0);
% 
% PsiF = V*Fi*PsiF;
% PhiF = V*Fi*PhiF;

% forcing mode
for k = 1:4
    figure(6+(k-1))
    contourf(xgrid',ygrid',reshape(real(PhiF(:,k)),[nx,ny]), 64,'LineStyle','none')
    %     colorbar; colormap(b2r(-0.02, 0.02));
    axis equal
%     title(sprintf('forcing mode %d, f %02f', k, omg))
    set(gca,'Fontname',font_name,'Fontsize',8,'LineWidth',0.5)
%     xlabel('$x$','Interpreter','Latex','fontsize',font_axes)
%     ylabel('$y$','Interpreter','Latex','fontsize',font_axes)
    plot_setting
     %     saveas(gcf, sprintf('./figures/forcing_%d_f_%02f.png', k, omg))
    print(gcf, sprintf('./figures/forcing_%d_f_%02f.png', k, omg), '-dpng')
end

% response mode
for k =1:4
    figure(10+(k-1))
    contourf(xgrid',ygrid',reshape(real(PsiF(:,k)),[nx,ny]), 64,'LineStyle','none')
%     colorbar; colormap(b2r(-0.02, 0.02));
    axis equal
%     title(sprintf('response mode %d, f %f', k, omg))
    set(gca,'Fontname',font_name,'Fontsize',8,'LineWidth',0.5)
%     xlabel('$x$','Interpreter','Latex','fontsize',font_axes)
%     ylabel('$y$','Interpreter','Latex','fontsize',font_axes)
    plot_setting
%     saveas(gcf, sprintf('./figures/response_%d_f_%02f.png', k, omg))
    print(gcf, sprintf('./figures/response_%d_f_%02f.png', k, omg), '-dpng')

end
