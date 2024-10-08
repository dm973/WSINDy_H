%% Z pinch data
ep = 0.15;%[0.1 0.27 0.6];
r2 = 1;
% load(['~/Desktop/learnJsigma_ep',num2str(ep),'_r2',num2str(r2),'.mat'],'traj','t','N')
load(['~/Desktop/Jsigma_data_ep',num2str(ep),'_r2',num2str(r2),'.mat'],'traj','t','N')
N = length(traj)
tol = mean(diff(t));
% load('~/Desktop/J_sigma_ep01_singletraj_results.mat','Ws');
% P =  dist(Ws~=0);
% [~,I] = sort(P(:,1),'ascend');
% inds = I(1:end);
inds = 1:N;

Ws = []; errs_J = []; errs_J0 = [];
tic;
Nmax = length(inds);
X_sigma = cell(Nmax,1);
t_sigma = cell(Nmax,1);
tcap = 20000;
for nn=1:Nmax
    [tnew,y0,iout,jout] = intersections(t(1:tcap),t(1:tcap)*0,t(1:tcap),traj{inds(nn)}(1:tcap,4),0);
    disp([nn toc])
    Xtemp = interp1(t(1:tcap),traj{inds(nn)}(1:tcap,1:3),tnew);    
    yy = find(Xtemp(:,2)>=0);
    X_sigma{nn} = [real(Xtemp(yy,[1 3])) sqrt(r2)*ones(length(yy),1)];
    t_sigma{nn} = tnew(yy)*ep;
end

%%%
B0 = 1; a1 = 0.3; a2 = 0.3; kx1 = 3; ky1 = 1; kx2 = 1; ky2 = 3;
B = @(q1,p1,q2,p2) B0 + a1*cos(kx1*q1+ky1*q2)+a2*cos(kx2*q1+ky2*q2);
B = @(x) B(x(1),0,x(2),0)/ep;
Bvar = @(x,y,r)B([x y]);
Bgrad_x = @(x) -a1*kx1*sin(kx1*x(1)+ky1*x(2))-a2*kx2*sin(kx2*x(1)+ky2*x(2))/ep;
Bgrad_y = @(x) -a1*ky1*sin(kx1*x(1)+ky1*x(2))-a2*ky2*sin(kx2*x(1)+ky2*x(2))/ep;
Bgrad = @(x) [-a1*kx1*sin(kx1*x(1)+ky1*x(2))-a2*kx2*sin(kx2*x(1)+ky2*x(2)),... 
    -a1*ky1*sin(kx1*x(1)+ky1*x(2))-a2*ky2*sin(kx2*x(1)+ky2*x(2)), 0]/ep;

% term('fHandle',@(x,y,r)B([x y]).^-1)
% Bgradterms = [term('fHandle',@(x1,x2,x3)Bgrad_x(x1 x2))]

dfuns = {@(x)0,@(x)1,@(x)B(x)/x(3)};
Afuns = {@(x)0,@(x)1,@(x)0;@(x)-1,@(x)0,@(x)0;@(x)0,@(x)0,@(x)0};

%% build lib 
addpath(genpath('../wsindy_obj_base'))

% lib = library('nstates',3);
% lib.add_terms(term('fHandle',@(x,y,r)r.^2/B([x y])/2));
% % lib.add_terms(term('fHandle',@(x,y,r)(r.^3/2)/B([x y]).^2.*Bgrad_y([x y])));
% % lib2 = get_trigprod_lib(3,[1 2 3 4 6 9],1,[],0,0);
% lib2 = get_trigprod_lib(2,[1 2 3 4 6 8 9 12],1,[],1,0);
% tags = cellfun(@(t)[t 1],lib2.tags,'un',0);
% lib.add_terms(tags);
% % tags = get_tags(1:3,[],3);
% % tags = tags(tags(:,3)<=1,:);
% % tags = prodtags(tags,[]);
% % lib.add_terms(tags);

lib = library('nstates',3);
lib2 = get_trigprod_lib(2,[0 kx1:kx1:3*kx1 ky1:ky1:3*ky1],1,[],1,0);
tags = cellfun(@(t)[t 2],lib2.tags,'un',0);
Bpows= [-1 -3];
Bpowterms = arrayfun(@(i)getBpowterm(Bvar,i,3),Bpows);
for j=1:length(Bpows)
    cellfun(@(tt)lib.add_terms(prodterm(Bpowterms(j),term('ftag',tt,'coeff',0.5))),tags);
end

lib = library('nstates',3);
lib2 = library('nstates',3);
fmax = 5;
freq = [[1 3];[3 1]];
tags = trigtags(freq,fmax);
tags = prodtags(tags,[]);
cellfun(@(tt)lib2.add_terms(term('ftag',[tt 2],'coeff',0.5,'gradon',1)),tags,'un',0);
Bpows= [-1 -3 -5];
Bpowterms = arrayfun(@(i)getBpowterm(Bvar,i,3),Bpows);
for j=1:length(Bpows)
    lib.add_terms(cellfun(@(tt)prodterm(Bpowterms(j),tt,'gradon',1),lib2.terms,'un',0));
end

%% get G b

Jfuns = cellfun(@(tt) @(x)tt.evalterm(x), lib.terms, 'Un',0);
Jgradfuns = cellfun(@(tt) @(x)cell2mat(tt.evalgrads(x)), lib.terms, 'Un',0);

%%% build G
G = []; Jmat = []; D = [];
subsamp = 4;
for ind=1:Nmax
    [xdot,x,~,~] = getxdot(X_sigma{ind},t_sigma{ind});
    x = x(1:subsamp:end,:);
    M = size(x,1);
    xdot = xdot(1:subsamp:end,:);
    [G_temp,Jcell,denom,Acell,P] = buildG(x,xdot,Jgradfuns,dfuns,Afuns);
    D_temp = cell2mat(cellfun(@(fn)dot(fn',denom,2),Jcell,'un',0)');
    G = [G;G_temp];
    D = [D;D_temp];
    Jmat = [Jmat;cell2mat(cellfun(@(f)arrayfun(@(i)f(x(i,:)),(1:M)'),Jfuns,'un',0))];
    disp([ind toc])
end

% G = [[J0] ----- [xdot*(d*gradJ_i) - A*gradJ_i] -----]
% min||G*w||^2 s.t. w(1) = 1, D*w > 0
% min||G*w||^2 s.t. ||w|| = 1, D*w > 0 

%% enforcing J0 - stabilized way
tol = min(D(:,1))/2;
b = G(:,1);
A = G(:,2:end);
% S = spdiags([ones(size(G,1),3)],-2:0,size(G,1),size(G,1));
% S = S(:,1:3:end)*D(:,1);
% A = A./S;
% b = b./S;

Wsub = quadprog(A'*A,b'*A);
lambda = ep^3;
biginds = true(length(Wsub),1);
W = [1;Wsub];
%clf;plot(D*W)
%drawnow
maxits = length(Wsub);
max_constrow = 1000;

% ||w_iJ_i||/||J0|| >= lambda

for n=1:maxits
    % biginds_new = min(1,vecnorm(Jmat(:,2:end))/norm(Jmat(:,1))).*abs(Wsub')>=lambda;
    biginds_new = vecnorm(Jmat(:,2:end))/norm(Jmat(:,1)).*abs(Wsub')>=lambda;
    if ~isequal(biginds,biginds_new)
        biginds = biginds_new;
        % subD = true(size(D,1),1);
        % subD(randperm(length(subD),ceil((1-subfrac)*length(subD))))=0;
        [~,subD] = sort(D*W,'ascend');
        subD = subD(1:min(end,max_constrow));
        % inds = find(subD);
        % subD(inds(1:3:end)) = 0;
        % subD = subD(randperm(length(subD),ceil(subfrac*length(subD))));
        Wsub(~biginds) = 0;
        if any(Wsub)
            A_temp = A(:,biginds);
            Aineq_temp = -D(subD,1+find(biginds));
            bineq_temp = -tol*ones(length(subD),1)+D(subD,1);
            Wsub(biginds) = quadprog(A_temp'*A_temp,b'*A_temp,Aineq_temp,bineq_temp);
        end
    end
    W = [1;Wsub];
    %clf;plot(D*W)
    %drawnow
end

Ws = [Ws W];

%% get vector field
J_tol = 0;
r = sqrt(r2);
rhs_learned = @(x) rhs_J(x,Jgradfuns,W,dfuns,Afuns,J_tol);
rhs_learned_xy = @(x) nout(rhs_learned([x(1);x(2);r]),[1;2]);
rhs_learned_0 = @(x) rhs_J(x,Jgradfuns,[1;Wsub*0],dfuns,Afuns,J_tol);
rhs_learned_xy_0 = @(x) nout(rhs_learned_0([x(1);x(2);r]),[1;2]);

%% view vector field and time derivative
for ind=1:Nmax
    t_temp = t_sigma{ind};
    M = length(t_temp);
    [xdot,x,M,t_temp] = getxdot(X_sigma{ind},t_temp);
    
    V = arrayfun(@(x,y)rhs_learned_xy([x y]),x(:,1),x(:,2),'Un',0);
    V0 = arrayfun(@(x,y)rhs_learned_xy_0([x y]),x(:,1),x(:,2),'Un',0);
    err_x = [norm(xdot(:,1)-cellfun(@(v)v(1),V0)) norm(xdot(:,1)-cellfun(@(v)v(1),V))];
    err_y = [norm(xdot(:,2)-cellfun(@(v)v(2),V0)) norm(xdot(:,2)-cellfun(@(v)v(2),V))];
    
    
    errs_J0 = [errs_J0 [err_x(1);err_y(1)]];
    errs_J = [errs_J [err_x(2);err_y(2)]];
    
    disp([ind toc])
    disp([err_x err_y])
    
    figure(ind)
    subplot(2,2,[1 3])
    plot(xdot(:,1),xdot(:,2),cellfun(@(v)v(1),V),cellfun(@(v)v(2),V0), cellfun(@(v)v(1),V),cellfun(@(v)v(2),V))
    legend({'data','J0','learned'})
    
    subplot(2,2,2)
    plot(t_temp,xdot(:,1),t_temp,cellfun(@(v)v(1),V0), t_temp,cellfun(@(v)v(1),V))
    legend({'data','J0','learned'})
    title(err_x)
    
    subplot(2,2,4)
    plot(t_temp,xdot(:,2),t_temp,cellfun(@(v)v(2),V0), t_temp,cellfun(@(v)v(2),V))
    legend({'data','J0','learned'})
    title(err_y)
    drawnow
end

%% sim learned dynamics

x_Ls = cell(Nmax,1);
for ind = 1:Nmax
    t = t_sigma{ind};
    M = length(t);
    x = X_sigma{ind}(:,1:2);
    
    x0 = x(1,:);
    tol_dd = 10^-6;
    options_ode_sim = odeset('RelTol',tol_dd,'AbsTol',tol_dd*ones(length(x0),1));
    
    [t_L,x_Ls{ind}]=ode15s(@(t,x)rhs_learned_xy(x),t,x0,options_ode_sim);
    % [t_L,x_L]=ode15s(@(t,x)rhs_learned_xy_0(x),t,x0,options_ode_sim);
    
    figure(33);
    clf; 
    subplot(2,1,1)
    plot(t,x,t_L,x_Ls{ind},'--','linewidth',3)
    subplot(2,1,2)
    plot(x(:,1),x(:,2),x_Ls{ind}(:,1),x_Ls{ind}(:,2),'--','linewidth',3)
    drawnow
    disp([ind toc])
end

%% contour plot of J sigma and trajectories
addpath(genpath('../wsindy_obj_base'))
% J_sigma_ep01_multitraj_results.mat
% J_sigma_ep01_multitraj_results_trig.mat
% J_sigma_ep01_r2_all_7-20.mat
% J_sigma_ep01_r2_all_7-21.mat
% J_sigma_ep0.1_r2_all_7-22.mat
% J_sigma_ep01_r2_all.mat
% J_sigma_ep01_singletraj_results.mat
load('~/Desktop/J_sigma_ep01_r2_all.mat','lib','ep','Ws','X_sigma','t_sigma','x_Ls','Jfuns','r2')
Nmax = length(X_sigma);
if ~exist('r2','var')
    r2 = 2;
end
if ~exist('ep','var')
    ep = 0.1;
end
r = sqrt(r2);

W = Ws;
J = @(x) get_J(x,Jfuns,W);
Ng = 80;
ex = 4;
x_grid = linspace(min(cellfun(@(x)min(x(:,1))-ex,X_sigma(1:Nmax))),...
    max(cellfun(@(x)max(x(:,1))+ex,X_sigma(1:Nmax))),Ng);
y_grid = linspace(min(cellfun(@(x)min(x(:,2))-ex,X_sigma(1:Nmax))),...
    max(cellfun(@(x)max(x(:,2))+ex,X_sigma(1:Nmax))),Ng);
[xx,yy] = meshgrid(x_grid,y_grid);

capfun = @(x)atan(x/ep);
Q = arrayfun(@(x,y)capfun(J([x y r])),xx,yy);

ics = cell2mat(cellfun(@(x)x(1,1:2),X_sigma(:),'un',0));
nc = 40;

figure(Nmax+2);clf
contourf(xx,yy,Q,nc,'edgecolor','k')
hold on
for j=1:Nmax
    plot(X_sigma{j}(:,1),X_sigma{j}(:,2),'w','linewidth',3)
end
if exist('x_Ls','var')
    for j=1:Nmax
        plot(x_Ls{j}(:,1),x_Ls{j}(:,2),'r','linewidth',3)
    end
end
plot(ics(:,1),ics(:,2),'gx','markersize',6,'linewidth',3)
hold off
colorbar
% xlim([min(x_grid) max(x_grid)])
% ylim([min(y_grid) max(y_grid)])
L=get(gca,'children');
legend(L([end end-1 2]),{'J_\sigma','data','learned'},'Location','bestoutside');
set(gca,'fontsize',20)

%% save

% save(['~/Desktop/J_sigma_ep',num2str(ep),'_r',num2str(r2),'_all_7-22.mat'],'Ws','errs_J','errs_J0','lib','lambda','tol','X_sigma','t_sigma','x_Ls','B','Bgrad','dfuns','Afuns','Jfuns','Jgradfuns','max_constrow','J_tol')

% %% view phase plot vector field
% 
% Ng = 100;
% x_grid = linspace(min(cellfun(@(x)min(x(:,1)),X_sigma(1:Nmax))),...
%     max(cellfun(@(x)max(x(:,1)),X_sigma(1:Nmax))),Ng);
% y_grid = linspace(min(cellfun(@(x)min(x(:,2)),X_sigma(1:Nmax))),...
%     max(cellfun(@(x)max(x(:,2)),X_sigma(1:Nmax))),Ng);
% 
% % x_grid = linspace(1.75,3,Ng);
% % y_grid = linspace(2.25,2.75,Ng);
% 
% [xx,yy] = meshgrid(x_grid,y_grid);
% Q = arrayfun(@(x,y)rhs_learned_xy([x y]),xx,yy,'Un',0);
% 
% Qmax = inf;
% nc = 40;
% 
% figure(1)
% contourf(xx,yy,max(min(cellfun(@(q)q(1),Q),Qmax),-Qmax),nc)
% hold on
% for j=1:Nmax
%     plot(X_sigma{j}(:,1),X_sigma{j}(:,2),'w','linewidth',3)
% end
% if exist('x_Ls','var')
% for j=1:Nmax
%     plot(x_Ls{j}(:,1),x_Ls{j}(:,2),'r','linewidth',3)
% end
% end    
% ics = cell2mat(cellfun(@(x)x(1,1:2),X_sigma(:),'un',0));
% plot(ics(:,1),ics(:,2),'gx','markersize',6,'linewidth',3)
% hold off
% colorbar
% xlim([min(x_grid) max(x_grid)])
% ylim([min(y_grid) max(y_grid)])
% 
% figure(2)
% contourf(xx,yy,max(min(cellfun(@(q)q(2),Q),Qmax),-Qmax),nc)
% hold on
% for j=1:Nmax
%     plot(X_sigma{j}(:,1),X_sigma{j}(:,2),'w','linewidth',3)
% end
% if exist('x_Ls','var')
% for j=1:Nmax
%     plot(x_Ls{j}(:,1),x_Ls{j}(:,2),'r','linewidth',3)
% end
% end
% plot(ics(:,1),ics(:,2),'gx','markersize',6,'linewidth',3)
% hold off
% colorbar
% xlim([min(x_grid) max(x_grid)])
% ylim([min(y_grid) max(y_grid)])

%% get coeffs - naive way

% lambda = 0.001;
% W = sparsifyJsigma(G,lambda);
% W
% subsamp = 1;
% [U,S,V] = svd(G(randperm(size(G,1),ceil(size(G,1)*subsamp)),:),0);
% W = V(:,end);

%% get coeffs - stabilized way

% subD = find(D(:,1)<0.2*range(D(:,1)));
% subfrac = 0.1;
% subD = subD(randperm(length(subD),ceil(subfrac*length(subD))));
% tol = 10^-6;
% W = quadprog(G'*G,zeros(size(G,2),1),-D(subD,:),-tol*ones(length(subD),1));
% lambda = exp(mean(log(abs(W(2:end)))));
% biginds = true(length(W),1);
% excl_inds = 1;
% 
% for n=1:length(W)
%     biginds_new = abs(W)>=lambda;
%     biginds_new(excl_inds)=true;
%     if ~isequal(biginds,biginds_new)
%         biginds = biginds_new;
%         W(~biginds) = 0;
%         G_temp = G(:,biginds);
%         D_temp = D(subD,biginds);
%         W(biginds) = quadprog(G_temp'*G_temp,zeros(size(G_temp,2),1),-D_temp,zeros(length(subD),1));
%     end
% end

%% enforcing J0 - naive way
% b = G(:,1);
% A = G(:,2:end);
% lambda = 0.1;
% Wtemp = sparsifyDynamics(A,b,lambda,1,0,ones(size(A,2),1),inf,1,[]);
% W = [1;-Wtemp];

function t = getBpowterm(B,pow,nstates)
    Bterm = term('fHandle',B);
    Bgrads = repelem(term('gradon',0,'nstates',nstates),1,nstates);
    for j=1:nstates
        Bgrads(j) = term('nstates',nstates,'fHandle',@(varargin)pow*B(varargin{:}).^(pow-1).*Bterm.gradterms(j).fHandle(varargin{:}));
    end
    t = term('nstates',nstates,'fHandle',@(varargin)B(varargin{:}).^pow,'gradterms',Bgrads);
end

function W = sparsifyJsigma(G,lambda)
    check = 1;
    W = ones(size(G,2),1);
    while check
        S = find(W);
        % [~,~,VD] = svd(D(:,S),0);
        % [~,~,V] = svd(G(:,S)*VD,0);
        % c = V(:,end)
        % W(S) = VD*c
        [~,~,V] = svd(G(:,S),0);
        W(S) = V(:,end);
        smallinds = abs(W(S)) < lambda;
        if any(smallinds)
            W(S(smallinds)) = 0;
            check = 1;
        else
            check = 0;
        end
    end
end

function dx = rhs(x,w,a,c)
%%% J=1/2(x^2+y^2), A = [0 sin(w*x);-sin(w*x) 0], d = [x y];
   d = x(1)^2+x(2)^2;
   dx = x*0;
   dx(1) = (a+c*sin(w*x(1)))*x(2)/d;
   dx(2) = (-a-c*sin(w*x(1)))*x(1)/d;
end

function J = get_J(x,Jfuns,w)
    inds = find(w)';
    J = 0;
    for j = inds
        J = J + w(j)*Jfuns{j}(x);
    end
end

function dx = rhs_J(x,Jgradfuns,w,dfuns,Afuns,tol)
    d = length(x);
    Jvec = zeros(d,1);
    inds = find(w)';
    for j=inds
        dJ = w(j)*Jgradfuns{j}(x);
        Jvec = Jvec+dJ(:);
    end
    dvec = cellfun(@(f)f(x),dfuns(:));
    Amat = cellfun(@(f)f(x),Afuns);
    dx = (Amat*Jvec);
    dd = dot(dvec,Jvec);
    s = sign(dd);
    dd = s*max(abs(dd),tol);
    dx = dx/dd;
end

function [G,Jcell,denom,Acell,P] = buildG(x,xdot,Jgradfuns,dfuns,Afuns)
    Jcell = Jgradfun(x,Jgradfuns);
    denom = dfun(x,dfuns);
    Acell = Afun(x,Afuns);
    [M,d] = size(x);
    J = length(Jgradfuns);
    P = cell(M,1);
    for m=1:M
        P{m} = xdot(m,:)'*denom(m,:)-Acell{m};        
    end
    % D = cell2mat(cellfun(@(J)dot(denom,J',2),Jcell(:)','Un',0));
    P = blkdiag(P{:});
    G = zeros(M*d,J);
    for j=1:J
        G(:,j) = P*Jcell{j}(:);
    end
end

function Jcell = Jgradfun(x,Jgradfuns)
    [M,d] = size(x);
    J = length(Jgradfuns);
    Jcell = cell(J,1);
    for j=1:J
        Jmat = zeros(d,M);
        for m=1:M
            Jmat(:,m) = Jgradfuns{j}(x(m,:));
        end
        Jcell{j} = Jmat;
    end
end

function [xdot,x,M,t] = getxdot(x,t)
    [M,d] = size(x);
    xdot = zeros(M-2,d);
    for j=1:M-2
        xdot(j,:) = [(t(j+2) - t(j+1))*(x(j+1,:)-x(j,:))/(t(j+1) - t(j))+...
            (t(j+1) - t(j))*(x(j+2,:)-x(j+1,:))/(t(j+2) - t(j+1))]/(t(j+2) - t(j)); %%% fix
    end
    x = x(2:end-1,:);
    M = M-2;
    t = t(2:end-1);
end

function denom = dfun(x,dfuns)
    [M,d] = size(x);
    denom = zeros(M,d);
    for m=1:M
        for i=1:d
            denom(m,i) = dfuns{i}(x(m,:));
        end
    end
end

function Acell = Afun(x,Afuns)
    [M,d] = size(x);
    Acell = cell(M,1);
    for m=1:M
        A = zeros(d);
        for i=1:d
            for j=1:d
                A(i,j) = Afuns{i,j}(x(m,:));
            end
        end
        Acell{m} = A;
    end
end

function v=nout(v,inds)
    v = v(inds);
end

% tags = get_tags(1:2:5,[],nstates);
    % arrayfun(@(i) @(x)cos(tags(i,1)*x(1))*cos(tags(i,2)*x(2)),1:size(tags,1),'un',0 ),...
    % arrayfun(@(i) @(x)cos(tags(i,1)*x(1))*sin(tags(i,2)*x(2)),1:size(tags,1),'un',0 ),...
    % arrayfun(@(i) @(x)sin(tags(i,1)*x(1))*cos(tags(i,2)*x(2)),1:size(tags,1),'un',0 )];

    % Jfuns = {@(x)x(1).^2,@(x)x(2).^2,@(x)x(1).*x(2),@(x)x(3)^2/B([x(1) x(2)])/2};
% Jgradfuns = {@(x)[2*x(1) 0 0],@(x)[0 2*x(2) 0],@(x)[x(2) x(1) 0],...
%     @(x)-x(3)^2/2/B([x(1) x(2)])^2*Bgrad(x)+[0 0 x(3)/B([x(1) x(2)])]};
