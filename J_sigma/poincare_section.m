addpath(genpath('wsindy_obj_base'))
addpath(genpath('../'))

%%% E potential
f_E=[1 2 3];
coeff_E=1./(1:length(f_E));
f_E=[0];
coeff_E=0;
phi_E = arrayfun(@(f)@(q1,p1,q2,p2) sin(f*q1).*sin(f*q2),f_E,'uni',0);
phi_E_x = arrayfun(@(f)@(q1,p1,q2,p2) cos(f*q1).*sin(f*q2),f_E,'uni',0);
phi_E_y = arrayfun(@(f)@(q1,p1,q2,p2) sin(f*q1).*cos(f*q2),f_E,'uni',0);
phi_fun = @(q1,q2) q1*0;
for i=1:length(coeff_E)
    phi_fun = @(q1,q2) phi_fun(q1,q2) + coeff_E(i)*phi_E{i}(q1,0,q2,0);
end

%%% B field
% B0 = 1;
% B = @(q1,p1,q2,p2) B0+q1*0;
f_B=[1 2 3];
coeff_B=1./(1:length(f_B));
B_cell = arrayfun(@(f)@(q1,p1,q2,p2) sin(f*q1).*sin(f*q2),f_B,'uni',0);
B = @(q1,p1,q2,p2) q1*0;
for i=1:length(coeff_B)
    B = @(q1,p1,q2,p2) B(q1,p1,q2,p2) + coeff_B(i)*B_cell{i}(q1,0,q2,0);
end

B0 = 1; a1 = 0.3; a2 = 0.3; kx1 = 3; ky1 = 1; kx2 = 1; ky2 = 3;
B = @(q1,p1,q2,p2) B0 + a1*cos(kx1*q1+ky1*q2)+a2*cos(kx2*q1+ky2*q2);

%%% sim params
for r2 = [0.5 1 2 4]
for ep = [0.05 0.15 0.27]
Ny = 16; %%% num traj
% nc = 3000; np = 24;
% t = linspace(0,2*pi*nc,ceil(nc*np/2/pi)); %%% timespan
t = 0:0.001/ep^2:(50/ep^2);
tol_ode = 10^-10; %%% sim tol
X = [2 0] + linspace(0,2*pi,Ny)'*[0 1]; %%% diagonal X0 through pos space
X = [X;[0 0] + linspace(0,2*pi,Ny)'*[0 1]]; %%% diagonal X0 through pos space
X = [X;[-2 0] + linspace(0,2*pi,Ny)'*[0 1]]; %%% diagonal X0 through pos space
% X = rand(N,2)*pi;
phX = phi_fun(X(:,1),X(:,2)); %%% evaluate reduced hamiltonian
E = r2/2*ep^2; %%% fix energy
N = size(X,1);
V = repmat(sqrt(2*E/ep^2)*[1 0],N,1); %%% initialize on section

%%% dynamics
features = cell(4,1);
features{1} = {@(q1,p1,q2,p2) p1};
features{2} = [phi_E_x,{@(q1,p1,q2,p2)B(q1,p1,q2,p2).*p2}];
features{3} = {@(q1,p1,q2,p2) p2};
features{4} = [phi_E_y,{@(q1,p1,q2,p2)B(q1,p1,q2,p2).*p1}];
params = {[ep],[-f_E.*coeff_E 1],[ep],[-f_E.*coeff_E -1]};
rhs_p = @(x,params) rhs_fun(features,params,x);

traj = cell(N,1);
% Sigma_traj = cell(N,1);
parfor nn=1:N
    x0 = reshape([X(nn,:);V(nn,:)],[],1)';
    options_ode_sim = odeset('RelTol',tol_ode,'AbsTol',tol_ode*ones(1,length(x0)));
    [~,traj{nn}]=ode15s(@(t,x)rhs_p(x,params),t,x0,options_ode_sim);
    % Sigma_traj{nn} = [traj{nn}(:,[1 3])  hypot(traj{nn}(:,2),traj{nn}(:,4))];
    disp(nn)
end
% save(['~/Desktop/learnJsigma_ep',num2str(ep),'.mat'])
save(['~/Desktop/Jsigma_data_ep',num2str(ep),'_r2',num2str(r2),'.mat'])
end
end
%% 

% dr = '/home/danielmessenger/Dropbox/Boulder/research/data/WENDy_data/ode_data/'; 
% load([dr,'GC_ep0.04_vperp0.30803_Ephi0.18823.mat'])
% load([dr,'GC_ep0.02_vperp0.43563_Ephi0.18823.mat'])
ep = 0.6;%[0.1 0.27 0.6];
load(['~/Desktop/learnJsigma_ep',num2str(ep),'.mat'])
%vars = {'x_red','t','true_prod_tags','true_coeff','params','rhs_p','Tfs','V0','x_full'};

clf
X_sigma = cell(N,1);
X_sigma_interp = cell(N,1);
t_sigma = cell(N,1);
t_sigma_interp = cell(N,1);
tol=2*mean(diff(t));
dt_interp = length(t)/2000;

for nn=1:N
    tic;
    xx = find(abs(traj{nn}(:,4))<tol);    % [~,xx]=findpeaks(-abs(traj{nn}(:,4)));    
    Xtemp = traj{nn}(xx,1:3);    
    tnew = t(xx);

    % [tnew,y0,iout,jout] = intersections(t,t*0,t,traj{nn}(:,4));
    % Xtemp = interp1(t,traj{nn}(:,1:3),tnew);    
    yy = find(Xtemp(:,2)>=0);

    X_sigma{nn} = real(Xtemp(yy,[1 3]));
    t_sigma{nn} = tnew(yy);
    t_sigma_interp{nn} = t_sigma{nn}(1):dt_interp:t_sigma{nn}(end);

    X_sigma_interp{nn} = interp1(t_sigma{nn},X_sigma{nn},t_sigma_interp{nn});
    % PS{nn} = mod(PS{nn},2*pi);
    % scatter(X_sigma{nn}(:,1),X_sigma{nn}(:,2),'.')
    scatter(X_sigma_interp{nn}(:,1),X_sigma_interp{nn}(:,2),'.')
    hold on
    disp([nn toc])
end
plot(X(:,1),X(:,2),'ro','linewidth',3)
hold off


B0 = 1; a1 = 0.3; a2 = 0.3; kx1 = 3; ky1 = 1; kx2 = 1; ky2 = 3;
B = @(q1,p1,q2,p2) B0 + a1*cos(kx1*q1+ky1*q2)+a2*cos(kx2*q1+ky2*q2);
B =@(x) B(x(1),0,x(2),0);
Bgrad=@(x) [-a1*kx1*sin(kx1*x(1)+ky1*x(2))-a2*kx2*sin(kx2*x(1)+ky2*x(2)),... 
    -a1*ky1*sin(kx1*x(1)+ky1*x(2))-a2*ky2*sin(kx2*x(1)+ky2*x(2)), 0];

clf
for ind=1:N
    plot(X_sigma{ind}(:,1),X_sigma{ind}(:,2),'.')
    hold on
    % plot(PS_r{ind}(:,3),'o-')
    % plot(traj{ind}(:,1),'.')
    % plot(Sigma_traj{ind}(:,1),Sigma_traj{ind}(:,2),'.-')
    drawnow
end
% PS=cell2mat(PS);
% vperp = cell2mat(vperp);
% plot(vperp,'.')
% scatter(PS(:,1),PS(:,2),'.')
% save([dr,['GC_diag_ep',num2str(ep),'_E',num2str(E),'.mat']])
% saveas(gcf,['~/Desktop/poincare_sec_',num2str(ep),'.png'])