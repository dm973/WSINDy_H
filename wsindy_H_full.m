% data params
rng(1)
ep = 0.01;
expnum = 1;
noise_ratio = 0.01;
train_inds = 18;
ncyc = 4;
subsamp = 1;

% tf params
eta = 12;
g = @(x) 1./x;
phifun = @(x) exp(-eta*g(1-x.^2));
% phifun = optTFcos(1,0);
% phifun = optTFpoly(3,0);
mtmax_fac = 1;
tf_meth = 'direct';
toggle_calibrate_tf = 0;
subind = -6;
Kmax = inf;
TphiTf = 4;

% optim params
L0_meth = 'MSTLS';
alpha = 0.01;
lambdas = 10.^linspace(-3,-1,200);
reg_0_param = 10^-inf;
toggle_jointthresh = 0;
toggle_wendy = 0;
discrep_type = 1;
discrep_its = 0;
discrep_step = 0.75;
discrep_meth = 'MSTLS';

% plotting params
test_ratio = 0.5;
range_fac = 0;
tol_dd = 10^-8;
num_x0_its = 40;
tfinal_x0_frac = 1/8;
toggle_sim = 1;
pdeg_x0 = 0;
toggle_plot_loss=1;

dr = '/home/danielmessenger/Dropbox/Boulder/research/data/WENDy_data/ode_data/';
vars = {'x_full','t','Tfs','params_full','rhs_p_full','true_coeff','true_prod_tags','rhs_p','params'};

if expnum==1
    ode_name = ['coupledosc_ep',strrep(num2str(ep),'0.',''),'_mu=1_p1-p2-03-07'];
    load([dr,ode_name],vars{:});
    true_prod_tags_full = {[2 0 0 0],[0 2 0 0],[0 0 2 0],[0 0 0 2],[[-2i 0 1 0];[1 0 2i 0]],[[2i 0 1 0];[1 0 -2i 0]]};
    true_prod_tags_full = [true_prod_tags_full,{[zeros(2) true_prod_tags{end}]}];
    true_coeff_full = [1/2 1/2 ep/2 ep/2 ep ep]';
    true_coeff_full = [true_coeff_full;0]; 
    Ji = [[0 1 0 0];[-1 0 0 0];[0 0 0 1];[0 0 -1 0]];
elseif expnum==3
    ode_name=['GC_ep',num2str(ep),'_E1.001_8-17.mat'];
    load([dr,ode_name],vars{:});
    true_prod_tags_full = {[0 2 0 0],[0 0 0 2],[-1i 0 -1i 0],[-2i 0 -2i 0],[-3i 0 -3i 0]};
    true_coeff_full = [ep^2/2 ep^2/2 ep ep/2 ep/3]';
    Ji = inv([[0 ep -1 0];[-ep 0 0 0];[1 0 0 ep];[0 0 -ep 0]]).';
elseif expnum==4
    ode_name=['GC_2p_r_ep',num2str(ep),'_8-18_.mat'];
    load([dr,ode_name],vars{:});
    true_prod_tags_full = {[0 2 0 0],[0 0 0 2],[-1i 0 -1i 0],[-2i 0 -2i 0],[-3i 0 -3i 0]};
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
%%
tags = true_prod_tags_full;
lib = library('tags',tags);
Hlib = cellfun(@(tm) tm.fHandle, lib.terms, 'uni',0);
J = length(lib.terms);
[w_true,rhs_H_true,Hfun_true] = get_true_dyn_script(lib.tags,true_prod_tags_full,true_coeff_full,params_full,rhs_p_full);

tic; warning('off','MATLAB:dispatcher:UnresolvedFunctionHandle');warning('off','MATLAB:legend:IgnoringExtraEntries');warning('off','MATLAB:load:variableNotFound');warning('off','MATLAB:rankDeficientMatrix')

%%% format data
rng_seed = rng().Seed; rng(rng_seed);
Uobj = wsindy_data(x,t(:)).trimend(T_final).coarsen(subsamp).addnoise(noise_ratio,'seed',rng_seed);
dims = Uobj.dims;
    
%%% get tf
tf_param = floor((TphiTf*ppTfast/subsamp-1)/2);
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
    reg0 = rank(WS.G{1},norm(WS.G{1})*reg_0_param);
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

w_plot = WS.weights;
rhs_H_learned = get_rhs_H(Hlib,w_plot,nstates,WS.Ji);

tpr = tpscore(w_plot,w_true);
b_plot = WS.b{1};
G_plot = WS.G{1};
disp(['---------------------'])
disp(['M=',num2str(Uobj.dims)])
disp(['TPR=',num2str(tpr),'; sparsity=',num2str(length(find(w_plot))),'/',num2str(size(G_plot,2))])
disp(['E2=',num2str(norm(WS.weights-w_true)/norm(w_true))])
disp(['size(G)=',num2str(size(G_plot))])
disp(['cond(G)=10^',num2str(log10(cond(G_plot)))])

%% plot results
warning('off','MATLAB:legend:IgnoringExtraEntries')

Ttest = floor(min(T_final*test_ratio));
Utest = wsindy_data(x,t(:)).trimend(Ttest).coarsen(subsamp);


options_ode_sim = odeset('RelTol',tol_dd,'AbsTol',tol_dd*ones(1,nstates),'Events',@(T,Y)myEvent(T,Y,100*max(vecnorm([Utest.Uobs{:}],2,2))));

dt = mean(diff(Utest.grid{1}));
v = cellfun(@(v)v.Cfs{1}(1,:),WS.tf{1},'uni',0);
mt = cellfun(@(v)(length(v)-1)/2,v(:)');    
disp(['-------------',num2str(train_inds),'-----------'])
disp(['Tphi=',num2str(2*mt+1)])
disp(['Tfast / dt=',num2str(Tfast/dt)])
disp(['2mt / Tfast =',num2str((2*mt+1)*dt/Tfast)])
disp(['Tslow / 2mt=',num2str(Tslow./(2*mt+1))])
disp(['Tfinal / Tslow=',num2str(Uobj.dims/Tslow)])

x0_opt = {'meth','ind','ind',1};
x0_reduced = Utest.get_x0(rhs_H_learned,x0_opt{:});
[t_learned,xH0_learned]=ode15s(@(t,x)rhs_H_learned(x),Utest.grid{1},x0_reduced,options_ode_sim);
tplot = Utest.grid{1};
xobs_plot = cell2mat(Utest.Uobs);
xtrue_plot = xobs_plot;

%%

figure(1)
clf
for i=1:nstates
    subplot(nstates/2,2,i)
    plot(tplot,xobs_plot(:,i),t_learned, xH0_learned(:,i))
    hold on
    plot(Uobj.grid{1}(1:tf{1}{1}.rads*2+1),tf{1}{1}.Cfs{1}(1,:)/max(tf{1}{1}.Cfs{1}(1,:))*max(abs(xH0_learned(:,i))))
    title(num2str(norm(diff([xobs_plot(1:length(t_learned),i) xH0_learned(:,i)],[],2))/norm(xobs_plot(1:length(t_learned),i))))
    legend
    xlim([0 Uobj.grid{1}(5*(tf{1}{1}.rads*2+1))])
end

if toggle_plot_loss==1
    figure(2)
    semilogx(loss_wsindy(2,:),loss_wsindy(1,:),'o-')
    title(['lambda=',num2str(loss_wsindy(2,find(loss_wsindy(1,:)==min(loss_wsindy(1,:)),1)))])
    legend('Loss')
end

%%

% clz = cell(16,1);
% clz{1} = 'c-';
% clz{4} = 'y--';
% clz{16} = 'r-';
% 
% if TphiTf==1
%     foo = 1;
% else
%     foo = 0;
% end
% 
% if TphiTf==4
%     lw=4;
% else
%     lw=8;
% end
% T = find(Uobj.grid{1}>Tslow,1);
% figure(10)
% inds = 2:4;
% if foo==1
%     clf
%     plot3(x(1:T,inds(1)),x(1:T,inds(2)),x(1:T,inds(3)),'k','linewidth',2)
% end
% hold on
% plot3(xH0_learned(1:T,inds(1)),xH0_learned(1:T,inds(2)),xH0_learned(1:T,inds(3)),[clz{TphiTf}],'linewidth',lw)
% view([70.3740   25.1834])
% xlabel('$P$','interpreter','latex'), ylabel('$q$','interpreter','latex'), zlabel('$p$','interpreter','latex')
% 
% figure(11)
% inds = 1:3;
% if foo==1
%     clf
%     plot3(x(1:T,inds(1)),x(1:T,inds(2)),x(1:T,inds(3)),'k','linewidth',2)
% end
% hold on
% plot3(xH0_learned(1:T,inds(1)),xH0_learned(1:T,inds(2)),xH0_learned(1:T,inds(3)),[clz{TphiTf}],'linewidth',lw)
% view([100.8978   56.6380])
% xlabel('$Q$','interpreter','latex'), ylabel('$P$','interpreter','latex'), zlabel('$q$','interpreter','latex')
% 
% 
% clz = cell(16,1);
% clz{1} = 'c-';
% clz{4} = 'y-';
% clz{16} = 'r--';
% 
% lws = [2 2 3 3];
% legs = {'$Q$','$P$','$q$','$p$'};
% for i=1:nstates
%     figure(i);set(gcf,'position',[ 1756         393        1294         354])
%     if foo == 1
%         plot(Uobj.grid{1},Uobj.Uobs{i},'ko','linewidth',1.5)
%     end
%     hold on
%     if and(TphiTf==16,i<3)
%         clz{16} = 'r-';
%     else
%         clz{16} = 'r--';
%     end
%     plot(t_learned,xH0_learned(:,i),[clz{TphiTf}],'linewidth',lws(i))
%     xlim([0 Tslow])
% %     xlim([500 530])
%     ylim([-2.5 2.5])
%     if mod(i,2)==0
%         xlabel('$t$','interpreter','latex')
%     end
%     ylabel(legs{i},'interpreter','latex')    
%     set(gca,'fontsize',20,'ticklabelinterpreter','latex','position',[0.130000000000000   0.279625852909218   0.775000000000000   0.645374147090782])
%     grid on
%     if TphiTf==16
%         if i==3
%             L=get(gca,'children');
%             legend([L(4);L(3);L(2);L(1)],{'\bf Z','$\hat{\bf Z}; \sigma_{\phi f}=1$','$\hat{\bf Z};\sigma_{\phi f}=4$','$\hat{\bf Z};\sigma_{\phi f}=16$'},'interpreter','latex','fontsize',18,'location','best')
%         end
%         saveas(gcf,['~/Desktop/',num2str(i),'_full.png'])
%     end
% end
% 
% if TphiTf==16
%     for i=10:11
%         figure(i)
%         set(gca,'ticklabelinterpreter','latex','fontsize',14)
%         grid on
%         axis equal
%         if i==10
%             saveas(gcf,['~/Desktop/exp1_ind20_full_Pqp.png'])
%             L=get(gca,'children');
%             legend([L(4);L(3);L(2);L(1)],{'\bf Z','$\hat{\bf Z}; \sigma_{\phi f}=1$','$\hat{\bf Z};\sigma_{\phi f}=4$','$\hat{\bf Z};\sigma_{\phi f}=16$'},'interpreter','latex')
%         else
%             saveas(gcf,['~/Desktop/exp1_ind20_full_QPq.png'])
%         end
%     end
% end