disp(['-----loading data-----'])

dr = '/home/danielmessenger/Dropbox/Boulder/research/data/WENDy_data/ode_data/';
vars = {'x_red','t','true_prod_tags','true_coeff','params','rhs_p','Tfs','V0','x_full'};
if expnum==1
    ode_name =['coupledosc_ep',strrep(num2str(ep),'0.',''),'_mu=1.mat']; 
%     ode_name = ['coupledosc_ep',strrep(num2str(ep),'0.',''),'_mu=0.3.mat'];
    ode_name = ['coupledosc_ep',strrep(num2str(ep),'0.',''),'_mu=1_p1-p2-03-07'];
%     libtype = 'polytrigprod';
    plot_style = [];
elseif expnum==2
    Q0 = (2^mQ-1)/2^mQ*pi;
%     ode_name = ['pendulum_ep',strrep(num2str(ep),'0.',''),'_P=0-Q=pi',num2str(Q0/pi),'.mat'];
%     ode_name = ['pendulum_2ep',strrep(num2str(ep),'0.',''),'_P=0-Q=pi',num2str(Q0/pi),'.mat'];
    ode_name = ['pendulum_4ep',strrep(num2str(ep),'0.',''),'_P=0-Q=pi2',num2str(Q0/pi),'.mat'];
%     libtype = 'poly';
    plot_style = 'Henon-Heiles';
elseif expnum==3
    E = 1.346617691642832*ep;
%     ode_name = ['GC_ep',num2str(ep),'_E',num2str(E),'.mat'];
    ode_name = ['GC_ep',num2str(ep),'_E1.3465.mat'];
    ode_name = ['GC_ep',num2str(ep),'_E0.62356.mat'];
    ode_name=['GC_ep',num2str(ep),'_E1.001_8-17.mat'];

%     libtype = 'trigprod';
    plot_style = [];
elseif expnum==4
    q1 = pi; q2 = 0;
%     ode_name=['GC_2p_r_ep',strrep(num2str(ep),'0.',''),'_q1=',num2str(q1),'_q2=',num2str(q2),'.mat'];
%     libtype = 'trigprod';
    plot_style = [];
    ode_name = ['GC_2p_r_ep',num2str(ep),'_8-17_'];
    ode_name = ['GC_2p_r_ep',num2str(ep),'_8-18_'];
elseif expnum==5
    ode_name=['simple_osc_ep',num2str(ep),'_muep0',num2str(mu_ep),'.mat'];
%     libtype = 'polytrigprod';
    Tfrac=1;
    plot_style = [];
end
warning('off','MATLAB:load:variableNotFound')
load([dr,ode_name],vars{:});

% svd_learningorthogonaltransformation_nearlyperiodicsys;

%     [~,ii] = sort(Tfs(:,2)./Tfs(:,1),'descend');
ii = 1:length(x_red);
try 
    Tfast = Tfs(ii(test_inds),1);
    Tslow = Tfs(ii(test_inds),2);
catch
    Tfast = ii(test_inds)*0+2*pi;
    Tslow = Tfast/ep;
end

if and(ncyc>0,ppTf>0) %% set ppTf and ncyc directly
    dtnew = Tfast/ppTf;
    subsamp = ceil(dtnew/mean(diff(t)));
    T_final = min(floor(length(t)*Tslow*ncyc/range(t)),length(t));
elseif and(ncyc<0,ppTf>0) %% set ppTf directly, and -ncyc= desired number of total points
    dtnew = Tfast/ppTf;
    subsamp = ceil(dtnew/mean(diff(t)));
    T_final = ceil(-ncyc*subsamp);
elseif and(ncyc>0,ppTf<0) %% set ncyc directly, and -ppTf=desired number of total points
    T_final = min(floor(length(t)*Tslow/range(t)*ncyc),length(t));
    subsamp = floor(T_final/(-ppTf));
elseif and(ncyc<0,ppTf<0) %% -ncyc=desired number of points, -ppTf=subsamp rate
    T_final = floor(length(t)/-ncyc);
    subsamp = -ppTf*(test_inds*0+1);
end
T_final = max(T_final,min_tp*subsamp);
Uobj = arrayfun(@(i)wsindy_data(x_red{i},t(:)),ii(min(train_inds,end)));
Uobj = arrayfun(@(i)Uobj(i).trimend(T_final(i)),(1:length(Uobj))');
Uobj = arrayfun(@(i)Uobj(i).coarsen(subsamp(i)),(1:length(Uobj))');

rng_seed = rng().Seed; rng(rng_seed);
Uobj = arrayfun(@(i)Uobj(i).addnoise(noise_ratio,'seed',rng_seed),(1:length(Uobj))');

Ttest = floor(min(T_final*test_ratio));
Utest = arrayfun(@(i)wsindy_data(x_red{i},t(:)),ii(min(test_inds,end)));
Utest = arrayfun(@(i)Utest(i).trimend(Ttest),(1:length(Utest))');
Utest = arrayfun(@(i)Utest(i).coarsen(max(subsamp)),(1:length(Utest))');

ntraj = length(Uobj);
nstates = Uobj(1).nstates;
dims = arrayfun(@(i)Uobj(i).dims,1:ntraj);

disp(['nstates=',num2str(nstates)])
disp(['num timepoints=',num2str(dims)])
disp(['train_inds=',num2str(train_inds(:)')])

%%% view data
if toggle_plot_sol~=0
    if toggle_plot_sol<0
        figure
    end
    subplot(3,1,1)
    Uobj(min(toggle_plot_sol,end)).plotDyn
    subplot(3,1,2)
    Uobj(min(toggle_plot_sol,end)).plotPhase
    drawnow
end