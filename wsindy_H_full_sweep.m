% data params
ep = 0.01;%
expnum = 1;%
noise_ratio = 0.01;
train_inds = 18;
ncyc = 4;%

% tf params
eta = 15;%
g = @(x) 1./x;%
phifun = @(x) exp(-eta*g(1-x.^2)+eta);%
mtmax_fac = 1;%
tf_meth = 'direct';%
toggle_calibrate_tf = 0;%
subind = -6;%
Kmax = inf;%

% optim params
L0_meth = 'MSTLS';%
alpha = 0.01;
reg_0_param = 10^-inf;
lambdas = 10.^linspace(-4,0,100);
toggle_jointthresh = 0;
toggle_wendy = 0;%
discrep_type = 1;%
discrep_its = 0;%
discrep_step = 0.75;%
discrep_meth = 'MSTLS';%

% plotting params
test_ratio = 0.25;%
range_fac = 0;%
tol_dd = 10^-8;%
num_x0_its = 40;%
tfinal_x0_frac = 1/8;%
toggle_sim = 1;%
pdeg_x0 = 0;%

dr = '/home/danielmessenger/Dropbox/Boulder/research/data/WENDy_data/ode_data/';
vars = {'x_full','t','Tfs','params_full','rhs_p_full','true_coeff','true_prod_tags'};

if expnum==1
    ode_name = ['coupledosc_ep',strrep(num2str(ep),'0.',''),'_mu=1_p1-p2-03-07'];
    load([dr,ode_name],vars{:});
    true_prod_tags_full = {[2 0 0 0],[0 2 0 0],[0 0 2 0],[0 0 0 2],[[-2i 0 1 0];[1 0 2i 0]],[[2i 0 1 0];[1 0 -2i 0]]};
    true_prod_tags_full = [true_prod_tags_full,{[zeros(2) true_prod_tags{end}]}];
    true_coeff_full = [1/2 1/2 ep/2 ep/2 ep ep]';
    true_coeff_full = [true_coeff_full;0];
    true_coeff_R0 = [1/2 1/2 0 0 0 0 0]';
    true_coeff_red = [0 0 ep/2 ep/2 0 0 true_coeff(end)]';
    Ji = [[0 1 0 0];[-1 0 0 0];[0 0 0 1];[0 0 -1 0]];
elseif expnum==3
    ode_name=['GC_ep',num2str(ep),'_E1.001_8-17.mat'];
    load([dr,ode_name],vars{:});
    true_prod_tags_full = {[0 2 0 0],[0 0 0 2],[-1i 0 -1i 0],[-2i 0 -2i 0],[-3i 0 -3i 0]};
%     true_prod_tags_full = [true_prod_tags_full,{[zeros(2) true_prod_tags{end}]}];
    true_coeff_full = [ep^2/2 ep^2/2 ep ep/2 ep/3]';
    Ji = inv([[0 ep -1 0];[-ep 0 0 0];[1 0 0 ep];[0 0 -ep 0]]).';
end

Tfast = Tfs(train_inds,1);
Tslow = Tfs(train_inds,2);
dt = mean(diff(t));
ppTfast = Tfast/dt;
ppTslow = Tslow/dt;
T_final = ceil(ncyc*ppTslow);

x = x_full{train_inds};
nstates = size(x,2);
ntraj = 1;

tags = true_prod_tags_full;
lib = library('tags',tags);
Hlib = cellfun(@(tm) tm.fHandle, lib.terms, 'uni',0);
J = length(lib.terms);
[w_true,rhs_H_true,Hfun_true] = get_true_dyn_script(lib.tags,true_prod_tags_full,true_coeff_full,params_full,rhs_p_full);

%%% sweep

subsamps = 1;%:2*floor(ppTfast);
TphiTfs = 2.^(-1:0.25:5);
rngs = 1:100;

% TphiTfs = TphiTfs(10:15);
% rngs = 1:24;

[mm,nn,zz] = ndgrid(subsamps,TphiTfs,rngs);
mm = mm(:);nn = nn(:);zz=zz(:);
errs = zeros(length(mm),1);
mods = cell(length(mm),1);

parfor kk=1:length(mm)
    tic; warning('off','MATLAB:dispatcher:UnresolvedFunctionHandle');warning('off','MATLAB:legend:IgnoringExtraEntries');warning('off','MATLAB:load:variableNotFound');warning('off','MATLAB:rankDeficientMatrix')

    %%% format data
    subsamp = mm(kk);
    rng(zz(kk)); rng_seed = rng().Seed;
    Uobj = wsindy_data(x,t(:)).trimend(T_final).coarsen(subsamp).addnoise(noise_ratio,'seed',rng_seed);
    dims = Uobj.dims;
    
    %%% get tf
    TphiTf = nn(kk);
    tf_param = floor((TphiTf*ppTfast/subsamp-1)/2)    
    if and(isequal(tf_meth,'direct'),tf_param<0)
        arrayfun(@(U)U.getcorners,Uobj);
        tf_param0 = get_tf_support(phifun,Uobj(1).dims,-tf_param,mean(arrayfun(@(j)1/mean(1./arrayfun(@(i)Uobj(j).ks(i),1:nstates)),1:ntraj)));
        tf = {arrayfun(@(j)testfcn(Uobj,'stateind',j,'meth','direct','param',tf_param0,'phifuns',phifun,'Kmax',Kmax,'mtmax',floor((dims-1)*mtmax_fac/2),'subinds',subind),(1:nstates)','uni',0)};
    else
        tf = {arrayfun(@(j)testfcn(Uobj,'stateind',j,'meth',tf_meth,'param',tf_param,'phifuns',phifun,'Kmax',Kmax,'mtmax',floor((dims-1)*mtmax_fac/2),'subinds',subind),(1:nstates)','uni',0)};
    end
    
    %%% solve linear system
    WS = wsindy_model(Uobj,lib,tf,'toggleH',1,'Ji',Ji);
    WS.cat_Gb;
    optm = WS_opt();
    if toggle_calibrate_tf~=0
        WS = optm.ols_tf(WS,'res_tol',toggle_calibrate_tf,'m_inc',0.25);
    end
    if isequal(L0_meth,'MSTLS')
        % reg0 = rank(WS.G{1},norm(WS.G{1})*reg_0_param);
        [WS,loss_wsindy,its] = optm.MSTLS(WS,'toggle_jointthresh',toggle_jointthresh,...
        'reg0', 0, 'alpha',alpha,'lambdas',lambdas);
    elseif isequal(L0_meth,'ols')
        WS = optm.ols(WS);    
    elseif isequal(L0_meth,'subspace_pursuit')
        s = rank(WS.G{1},norm(WS.G{1})*reg_0_param);
        ncv = min(floor(length(lib(1).terms)/2),100);
        WS = optm.subspacePursuitCV(WS,s,ncv,0);
    end
    if toggle_wendy > 0
        [WS,w_its,res,res_0,CovW] = optm.wendy(WS,'maxits',toggle_wendy);
    end
    for k=1:discrep_its
        for i=1:WS.ntraj
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
            s_new = s - length(find(WS.weights)); 
            ncv = floor(length(tags)/2);
            WS = optm.subspacePursuitCV(WS,s_new,ncv,discrep_type);
        end
    end
    mods{kk} = WS.weights;
    errs(kk) = max(abs(WS.weights-true_coeff_full)./abs(true_coeff_full));
    disp([kk toc])

end

%%%

errs_full = squeeze(mean(reshape(cellfun(@(W) norm(W-true_coeff_full)/norm(true_coeff_full),mods),...
    length(subsamps),length(TphiTfs),length(rngs)),3))';

errs_R0 = squeeze(mean(reshape(cellfun(@(W) norm(W-true_coeff_R0)/norm(true_coeff_R0),mods),...
    length(subsamps),length(TphiTfs),length(rngs)),3))';

errs_red = squeeze(mean(reshape(cellfun(@(W) norm(W-true_coeff_red)/norm(true_coeff_red),mods),...
    length(subsamps),length(TphiTfs),length(rngs)),3))';


tp_full = squeeze(mean(reshape(cellfun(@(W) tpscore(W,true_coeff_full),mods),...
    length(subsamps),length(TphiTfs),length(rngs)),3))';

tp_R0 = squeeze(mean(reshape(cellfun(@(W) tpscore(W,true_coeff_R0),mods),...
    length(subsamps),length(TphiTfs),length(rngs)),3))';

tp_red = squeeze(mean(reshape(cellfun(@(W) tpscore(W,true_coeff_red),mods),...
    length(subsamps),length(TphiTfs),length(rngs)),3))';

figure(11)
loglog(TphiTfs,min(1,[errs_R0 errs_full errs_red]),'o-','linewidth',3)
% xlim([1 16])
xlabel('$\sigma_{\phi f}$','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','FontSize',16,'Xtick',2.^(0:4))
legend({'$\Delta {\bf w}_0$','$\Delta {\bf w}_{full}$','$\Delta {\bf w}_{red}$'},...
    'Interpreter','latex','location','best')
grid on

figure(10)
semilogx(TphiTfs,[tp_R0 tp_full tp_red],'o-','linewidth',3)
% xlim([1 16])
xlabel('$\sigma_{\phi f}$','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','FontSize',16,'Xtick',2.^(0:4))
legend({'TPR$_0$','TPR$_{full}$','TPR$_{red}$'},'Interpreter','latex','location','best')
grid on

%%

errs2 = cellfun(@(W) max(abs(W-true_coeff_full)./abs(true_coeff_full)),mods);
% errs2 = cellfun(@(W) norm(W-true_coeff_full,inf)/norm(true_coeff_full),mods);

surf(reshape(ppTfast./mm,length(subsamps),[]),reshape(nn,length(subsamps),[]),...
    reshape(errs2,length(subsamps),[]));view([0 90])
colorbar
set(gca,'Xscale','log','Zscale','log')
xlabel('T_f/T_{\Delta t}')
ylabel('T_\phi/T_f')

%%

clf
lib = library('tags',true_prod_tags_full);
true_coeff_full(end) = true_coeff(end);
for ind = 1:length(true_coeff_full)
    subplot(4,2,ind)
    dat = cellfun(@(W)W(ind),reshape(mods,length(subsamps),[]));
    surf(reshape(ppTfast./mm,length(subsamps),[]),reshape(nn,length(subsamps),[]),...
        abs(dat-true_coeff_full(ind))/true_coeff_full(ind),'edgecolor','none');view([0 90])
    set(gca,'Xscale','log')
    colorbar
    xlabel('T_f/T_{\Delta t}')
    ylabel('T_\phi/T_f')
    title(functions(lib.terms{ind}.fHandle).function)
end

%%
clf
lib = library('tags',true_prod_tags_full);
% true_coeff_full(end) = true_coeff(end);
subs = [20];
Ji_r = [[0 0 1 0];[0 0 0 0];[-1 0 0 0];[0 0 0 0]];
for ind = 1:length(true_coeff_full)
    subplot(4,2,ind)
    dat = cellfun(@(W)W(ind),reshape(mods,length(subsamps),[]));
    dat = abs(dat(subs,:)'+true_coeff_full(ind)*ep)/abs(ep*true_coeff_full(ind));
    semilogy(TphiTfs,dat)
    xlabel('T_\phi/T_f')
    title(functions(lib.terms{ind}.fHandle).function)
    grid on
end
legend(arrayfun(@(x)num2str(x),ppTfast./subsamps(subs),'uni',0),'location','best')
