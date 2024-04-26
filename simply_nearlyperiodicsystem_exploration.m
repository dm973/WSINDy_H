tol_dd = 10^-10;
opt_4 = odeset('RelTol',tol_dd,'AbsTol',tol_dd*ones(1,4));
opt_2 = odeset('RelTol',tol_dd,'AbsTol',tol_dd*ones(1,2));

A = 5;
omeg = 2;
mu0 = @(Q,P,q,p) 0.5*(Q.^2+P.^2);
mu1 = @(Q,P,q,p) A/4*cos(omeg*q).*(P.^2-Q.^2);

for ep = 0.01:0.01:0.1
    t = linspace(0,6*2*pi/ep,10000);
    X0 = [1 0]; % fix to keep mu0 fixed
    q0 = 0; % fix to keep mu1 fixed
    R = 96;
    p0s = linspace(2,6,R);
    x0s = [ones(R,1)*[X0 q0] p0s'];
    mu_0=mu0(x0s(1,1),x0s(1,2),x0s(1,3),x0s(1,4));
    mu_1=mu1(x0s(1,1),x0s(1,2),x0s(1,3),x0s(1,4));
    mu_ep = mu_0+ep*mu_1;
    
    x_full = cell(R,1);
    x_red = cell(R,1);
    parfor r = 1:R
        [~,x_ep]=ode15s(@(t,x)rhs_ep(x,ep,A,omeg),t,x0s(r,:),opt_4);
        x_full{r} = x_ep;
        x_red{r} = x_ep(:,3:4);
        disp(r)
    end
    
    true_prod_tags = {[[2 0];[0 0]],[[0 2];[0 0]],[omeg*1i 0]};
    true_coeff = [ep/2 ep/2 -A*ep/2*mu_0];
    features = cell(2,1);
    features{1} = {@(q2,p2) p2};
    features{2} = {@(q2,p2) q2, @(q2,p2) sin(omeg*q2)};
    params = {ep,[-ep -omeg*A*ep/2*mu_0]};
    rhs_p = @(x,params) rhs_fun(features,params,x);
    Tfs = [2*pi*ones(R,1) 2*pi/ep*ones(R,1)];
    libtype = 'polytrigprod';
    dr = '/home/danielmessenger/Dropbox/Boulder/research/data/WENDy_data/ode_data/';
    save([dr,'simple_osc_ep',num2str(ep),'_muep0',num2str(mu_0),'.mat'])
end
%%

tol_dd = 10^-10;
opt_4 = odeset('RelTol',tol_dd,'AbsTol',tol_dd*ones(1,4));
opt_2 = odeset('RelTol',tol_dd,'AbsTol',tol_dd*ones(1,2));

A = 5;
omeg = 2;
mu0 = @(Q,P,q,p) 0.5*(Q.^2+P.^2);
mu1 = @(Q,P,q,p) A/4*cos(omeg*q).*(P.^2-Q.^2);

ep = 0.01;
t = linspace(0,6*2*pi/ep,10000);
X0 = [1 0]; % fix to keep mu0 fixed
q0 = 0; % fix to keep mu1 fixed
R = 96; 
p0s = linspace(2,6,R);
x0s = [ones(R,1)*[X0 q0] p0s'];
mu_0=mu0(x0s(1,1),x0s(1,2),x0s(1,3),x0s(1,4));
mu_1=mu1(x0s(1,1),x0s(1,2),x0s(1,3),x0s(1,4));
mu_ep = mu_0+ep*mu_1;

x0 = x0s(20,:);
[t_ep,x_ep]=ode15s(@(t,x)rhs_ep(x,ep,A,omeg),t,x0,opt_4);
[t_bar,x_bar]=ode15s(@(t,x)rhs_bar(x,ep,A,omeg),t,x0,opt_4);
[t_mu0,x_mu0]=ode15s(@(t,x)rhs_mu0(x,ep,mu_0,A,omeg),t,x0(3:4),opt_2);
[t_mu1,x_mu1]=ode15s(@(t,x)rhs_mu1(x,ep,mu_0,mu_1,A,omeg),t,x0(3:4),opt_2);

figure(1)
subplot(4,1,1)
plot(t_ep,x_ep(:,1),'k-',t_bar,x_bar(:,1),'g-*','markersize',1)
xlim([0 400])
legend({'coord1, full','coord1, Xav'})

subplot(4,1,2)
plot(t_ep,x_ep(:,2),'k-',t_bar,x_bar(:,2),'g-*','markersize',1)
xlim([0 400])
legend({'coord2, full','coord2, Xav'})

subplot(4,1,3)
plot(t_ep,x_ep(:,3),'k-',t_bar,x_bar(:,3),'g-*',t_mu0,x_mu0(:,1),'r--.',t_mu1,x_mu1(:,1),'c--','markersize',1)
title(['err 1st O=',num2str(rms(x_ep(:,3)-x_mu0(:,1))),'   err 2nd O=',num2str(rms(x_ep(:,3)-x_mu1(:,1)))])
legend({'coord3, full','coord3, Xav','coord3, H0red', 'coord3, H1red'})

subplot(4,1,4)
plot(t_ep,x_ep(:,4),'k-',t_bar,x_bar(:,4),'g-*',t_mu0,x_mu0(:,2),'r--.',t_mu1,x_mu1(:,2),'c--','markersize',1)
title(['err 1st O=',num2str(rms(x_ep(:,4)-x_mu0(:,2))),'   err 2nd O=',num2str(rms(x_ep(:,4)-x_mu1(:,2)))])
legend({'coord4, full','coord4, Xav','coord4, H0red', 'coord4, H1red'})
% 
% figure(2)
% plot(t_ep,mu0(x_ep(:,1),x_ep(:,2),x_ep(:,3),x_ep(:,4)),...
%     t_ep,mu0(x_ep(:,1),x_ep(:,2),x_ep(:,3),x_ep(:,4))+ep*mu1(x_ep(:,1),x_ep(:,2),x_ep(:,3),x_ep(:,4)))
% xlim([t_ep([1 end])])


function dz=rhs_ep(z,ep,A,omeg)
    dz = zeros(4,1);
    dz(1) = z(2);
    dz(2) = -z(1).*(1-A*ep*cos(omeg*z(3)));
    dz(3) = ep*z(4);
    dz(4) = ep*(-z(3)-omeg*A*0.5*z(1).^2.*sin(omeg*z(3)));
end

function dz=rhs_bar(z,ep,A,omeg)
    dz = zeros(4,1);
    dz(1) = z(2).*(1-A*ep/2*cos(omeg*z(3)));
    dz(2) = -z(1).*(1-A*ep/2*cos(omeg*z(3)));
    dz(3) = ep*z(4);
    dz(4) = ep*(-z(3)-omeg*A*1/4*(z(1).^2+z(2).^2).*sin(omeg*z(3)));
end

function dZ=rhs_mu0(Z,ep,mu0,A,omeg)
    dZ = zeros(2,1);
    dZ(1) = ep*Z(2);
    dZ(2) = ep*(-Z(1)-omeg*A*1/2*mu0.*sin(omeg*Z(1)));
end

function dZ=rhs_mu1(Z,ep,mu0,mu1,a,omeg)
    dZ = zeros(2,1);
    dZ(1) = ep*Z(2);
    dZ(2) = ep*(-Z(1)-omeg*a*1/2*(mu0+ep*mu1).*sin(omeg*Z(1)))...
        -omeg*(mu0+ep*mu1)/2*(ep*a*1/2)^2*sin(2*omeg*Z(1));
end
