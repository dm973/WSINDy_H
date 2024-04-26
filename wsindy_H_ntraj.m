addpath(genpath('wsindy_obj_base'))
%% set params
%%% load data params
clc;
rng('shuffle')

expnum = 3; % experiment number
ep = 0.05; % perturbation param
ppTf = 10; % points per fast mod6
ncyc = 4; % number of slow modes observed
min_tp = 25; % minimum number of timepoints
test_ratio = 1; % test on 1 slow
noise_ratio = 0; % noise ratio
train_inds = randperm(25,1)+5; % training inds
test_inds = [train_inds];
mu_ep = 0.5; % adiabatic invariant, expnum 5 only
mQ = 5; % pendulum proximity to unstable fp, expnum 2 only

use_presets = 1;

%%% lib params
toggle_get_lib = 1;
polys = [1:3];
trigs = [0:3]; L = 1;
include_crosstrig = 0;
add_polytrig = 1;
toggle_trimlib = 0; ord = 5; rangefun= @(b)b.^-2;
libtype = 'trigprod';

%%% tf params
eta = 9;
g = @(x) 1./x; % -log(x)
phifun = @(x) exp(-eta*g(1-x.^2))/exp(-eta); 
% phifun = optTFpoly(5,0); % 
% phifun = optTFcos(2,1);
mtmax_fac = 0.75;
sobol_inds = [];
tf_meth = 'FFT';
tf_param = 2;%floor(Uobj(1).dims/40);
toggle_calibrate_tf = 0;
subind = [];

%%% optim params
reg_0_param = 10^-inf;
toggle_jointthresh = 1;
L0_meth = 'MSTLS';
toggle_wendy = 0;
discrep_type = 2;
discrep_its = 0;
discrep_step = 0.75;
discrep_meth = 'MSTLS';

%%% plot params
toggle_plot_sol = 0;
toggle_plot_results = 1;
toggle_plot_sims = 1;
toggle_plot_loss = 1;
range_fac = 0;
ep_fac = 20*ep; 
tol_dd = 10^-6;
num_x0_its = 40;
tfinal_x0_frac = 1/8;
pdeg_x0 = 3;

if use_presets
if expnum==1
    mQ=1;
    ep = 0.05;
    train_inds = 14;
    polys = 1:2;
    trigs = 0:2:6; L = 1;
    add_polytrig = 1;
    pdeg_x0 = 0;
elseif expnum==2
    train_inds = 20;
    mQ=5;
    ep = 0.05;
    polys = 1:4;
    trigs = []; L = 1;
    add_polytrig = 0;
    pdeg_x0 = 0;
elseif expnum==3
    mQ=1;
    ep = 0.03;
    train_inds = 7;
    polys = [];
    trigs = 0:4; L = 1;
    add_polytrig = 0;
    pdeg_x0 = 3;
elseif expnum==4
    mQ=1;
    ep = 0.03;
    train_inds = 3;
    polys = 1:2;
    trigs = 0:3; L = 1;
    add_polytrig = 0;
    pdeg_x0 = 3;
end
test_inds = [train_inds];
end

%%% load data
load_data_script;
%% get lib

if toggle_get_lib ~= 0
    get_lib_script;
end

%% get initial test functions

tf_script;

%% initialize model

WS = wsindy_model(repmat(Uobj,length(sobol_inds)+1,1),lib,tf,'toggleH',1);

%% solve linear system
warning('off','MATLAB:rankDeficientMatrix')

optm = WS_opt();
WS.cat_Gb;

if toggle_calibrate_tf~=0
    WS = optm.ols_tf(WS,'res_tol',toggle_calibrate_tf,'m_inc',0.25);
end

if isequal(L0_meth,'MSTLS')
    reg0 = rank(WS.G{1},norm(WS.G{1})*reg_0_param);
    disp(['keeping initial largest ',num2str(reg0),' coefficients'])
    [WS,loss_wsindy,its] = optm.MSTLS(WS,'toggle_jointthresh',toggle_jointthresh, 'reg0', reg0,'Hreg',0);
elseif isequal(L0_meth,'subspace_pursuit')
%     s = 10;max(10,floor(length(lib(1).terms)/8));
    s = rank(WS.G{1},norm(WS.G{1})*reg_0_param);
    ncv = min(floor(length(lib(1).terms)/2),100);
    WS = optm.subspacePursuitCV(WS,s,ncv,0);
end

if toggle_wendy > 0
    [WS,w_its,res,res_0,CovW] = optm.wendy(WS,'maxits',toggle_wendy);
end

%%% discrep correction
for k=1:discrep_its
    for i=1:ntraj
        for j=1:nstates
            disp([WS.tf{i}{j}.rads ceil(WS.tf{i}{j}.rads*discrep_step)])
            WS.tf{i}{j}.param = ceil(WS.tf{i}{j}.rads*discrep_step);
            WS.tf{i}{j}.meth = 'direct';
            WS.tf{i}{j}.Cfs = cell(WS.ndims,1);
            WS.tf{i}{j}.Cfsfft = cell(WS.ndims,1);
            WS.tf{i}{j}.get_rads;
        end
    end
    WS.get_Gb;

    reg0 = rank(WS.G{1},norm(WS.G{1})*10^-6);
    if isequal(discrep_meth,'MSTLS')
        [WS,loss_wsindy,its] = optm.MSTLS(WS,'toggle_discrep',discrep_type,'toggle_jointthresh',toggle_jointthresh, 'reg0', reg0);
    else
        s = 20;
        s_new = s - length(find(WS.weights)); 
        ncv = floor(length(tags)/2);
        WS = optm.subspacePursuitCV(WS,s_new,ncv,discrep_type);
    end
end

Hfun = @(varargin) varargin{1}*0;
inds = find(WS.weights);
for i=1:length(inds)
    Hfun = @(varargin) Hfun(varargin{:}) + WS.weights(inds(i))*Hlib{inds(i)}(varargin{:});
end

numiter=0;
for iter=1:numiter
%%% iterative H correction
G = []; b = [];
inds_0 = WS.weights~=0;
for ind = 1:ntraj
    H0 = Hfun(Uobj(ind).Uobs{:});
    gradH0 = lib.grad2mat(Uobj(ind).Uobs,inds_0)*WS.weights(inds_0);
    gradH0 = reshape(gradH0,[],nstates);
    A = [];
    for i=1:nstates
        A = [A spdiags(gradH0(:,i),0,size(gradH0,1),size(gradH0,1))];
    end
    Th = A*WS.get_theta{ind}(:,~ismember(1:length(WS.weights),find(WS.weights~=0)));
    tf2 = testfcn(H0-mean(H0),'meth','FFT','param',4,'subinds',-6);
    G = [G;arrayfunvec(Th, @(c)tf2.test(c,0), 1)];
    b = [b;arrayfunvec(H0, @(c)tf2.test(c,1), 1)];
end
WS2 = wsindy_model(Uobj,[],[]);
WS2.G = {G}; WS2.b = {b};
optm2 = WS_opt();
WS2 = optm2.MSTLS(WS2,'toggle_jointthresh',toggle_jointthresh, 'reg0', reg0);
WS.weights(~inds_0) = WS2.weights;
inds = setdiff(find(WS.weights),find(inds_0));
for i=1:length(inds)
    Hfun = @(varargin) Hfun(varargin{:}) + WS.weights(inds(i))*Hlib{inds(i)}(varargin{:});
end
end

%% display results
warning('off','MATLAB:dispatcher:UnresolvedFunctionHandle');
[w_true,rhs_H_true,Hfun_true] = get_true_dyn_script(tags,true_prod_tags,true_coeff,params,rhs_p);

w_plot = WS.weights;
rhs_H_learned = get_rhs_H(Hlib,w_plot,nstates,WS.Ji);

tpr = tpscore(w_plot,w_true);
b_plot = WS.b{1};
G_plot = WS.G{1};
disp(['---------------------'])
disp(['TPR=',num2str(tpr),'; sparsity=',num2str(length(find(w_plot))),'/',num2str(size(G_plot,2))])
disp(['E2=',num2str(norm(WS.weights-w_true)/norm(w_true))])
disp(['size(G)=',num2str(size(G_plot))])
disp(['cond(G)=10^',num2str(log10(cond(G_plot)))])

%% plot results
warning('off','MATLAB:legend:IgnoringExtraEntries')

if toggle_plot_results>0

    for traj_ind = 1:length(test_inds)
        options_ode_sim = odeset('RelTol',tol_dd,'AbsTol',tol_dd*ones(1,nstates),'Events',@(T,Y)myEvent(T,Y,5*max(vecnorm([Utest(traj_ind).Uobs{:}],2,2))));

        dt = mean(diff(Utest(traj_ind).grid{1}));
        v = cellfun(@(v)v.Cfs{1}(1,:),WS.tf{min(traj_ind,end)},'uni',0);
        mt = cellfun(@(v)(length(v)-1)/2,v(:)');    
        disp(['-------------',num2str(test_inds(traj_ind)),'-----------'])
        disp(['radius=',num2str(mt)])
        disp(['Tfast / dt=',num2str(Tfast(traj_ind)/dt)])
        disp(['2mt / Tfast =',num2str((2*mt+1)*dt/Tfast(traj_ind))])
        disp(['Tslow / 2mt=',num2str(Tslow(traj_ind)./(2*mt+1))])
        disp(['Tfinal / Tslow=',num2str(Uobj(min(traj_ind,end)).dims/Tslow(traj_ind))])


        figure(traj_ind+1)
        if toggle_plot_sims == 1
            tic,
            if pdeg_x0>0
                tfinal_x0 = floor(Utest(traj_ind).dims*tfinal_x0_frac);
                x0_opt = {'meth','grid_search','wd_test',pdeg_x0+1:ceil(tfinal_x0/num_x0_its):tfinal_x0,'polydeg',pdeg_x0,'Tf',floor(Utest(traj_ind).dims(1)/4)};
                x0_reduced = Utest(traj_ind).get_x0(rhs_H_learned,x0_opt{:});
                x0_true = Utest(traj_ind).get_x0(rhs_H_true,x0_opt{:});
            elseif pdeg_x0==0
                x0_opt = {'meth','ind','ind',1};
                x0_reduced = Utest(traj_ind).get_x0(rhs_H_learned,x0_opt{:});
                x0_true = Utest(traj_ind).get_x0(rhs_H_true,x0_opt{:});
            else 
                % x0_reduced = mean(x_red{test_inds(traj_ind)}(1:find(t>=Tfast,1),:));
                % x0_true = x0_reduced;
                x0_opt = {'meth','trap'};
                x0_reduced = Utest(traj_ind).get_x0(rhs_H_learned,x0_opt{:});
                x0_true = Utest(traj_ind).get_x0(rhs_H_true,x0_opt{:});
            end
            [~,xH0_learned]=ode15s(@(t,x)rhs_H_learned(x),Utest(traj_ind).grid{1},x0_reduced,options_ode_sim);
            [~,xH0_true]=ode15s(@(t,x)rhs_H_true(x),Utest(traj_ind).grid{1},x0_true,options_ode_sim);
            disp(['time to sim=',num2str(toc)])
        end

        tplot = Utest(traj_ind).grid{1};
        xobs_plot = cell2mat(Utest(traj_ind).Uobs);
%         if ~isempty(Utest(traj_ind).noise{1})
%             xtrue_plot = cell2mat(Utest(traj_ind).Uobs) - cell2mat(Utest(traj_ind).noise(:,1)');
%         else
            xtrue_plot = xobs_plot;
%         end

        if ismember(test_inds(traj_ind),train_inds)
            v = cellfun(@(v)v.Cfs{1}(1,:),tf{min(find(test_inds(traj_ind)==train_inds),end)},'uni',0);
            mt = cellfun(@(v)(length(v)-1)/2,v(:)');
            kx = Uobj(test_inds(traj_ind)==train_inds).ks;
%             ks_vfft = [0:Uobj{ismember(test_inds(traj_ind),train_inds)}.dims-1]-ceil(Uobj{ismember(test_inds(traj_ind),train_inds)}.dims/2);
            vfft = abs(fftshift(cell2mat(arrayfun(@(i)fft([v{i}'/norm(v{i},1);zeros(length(tplot)-2*mt(i)-1,1)]),1:nstates,'uni',0))'));
        else
            kx = [];
            vfft = [];
        end
        if toggle_plot_sims == 1
            plot_wsindy_hamiltonian;
        else
            plot_no_sims;
        end
        drawnow
%         saveas(gcf,['~/Desktop/',ode_name,'_',num2str(traj_ind),'_MSMupdate_07-20-23_mtref.png'])
    end
end

if toggle_plot_loss==1
    figure(length(test_inds)+1)
    semilogx(loss_wsindy(2,:),loss_wsindy(1,:),'o-','linewidth',2)
    hold on
    plot( ...
        ...loss_wsindy(2,find(loss_wsindy(1,:)==min(loss_wsindy(1,:)),1)),min(loss_wsindy(1,:)),'kx',...
        loss_wsindy(2,find(loss_wsindy(1,:)==min(loss_wsindy(1,:)),1)),min(loss_wsindy(1,:)),'ko',...
        'markersize',20,'linewidth',4)
    hold off
    legend({'$\mathcal{L}(\lambda)$','$\widehat{\lambda}$'},'interpreter','latex')
    xlabel('$\lambda$','interpreter','latex')
    grid on
    set(gca,'fontsize',16,'ticklabelinterpreter','latex')
end

% figure(15)
% range(Hfun(x_red{train_inds(1)}(:,1),x_red{train_inds(1)}(:,2)))
% plot(Hfun(x_red{train_inds(1)}(:,1),x_red{train_inds(1)}(:,2)))
% i=20;

% semilogy(abs(fftshift(fft(Hfun_true(x_red{i}(:,1),x_red{i}(:,2))))))