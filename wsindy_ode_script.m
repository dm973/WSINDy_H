addpath(genpath('wsindy_obj_base'))
%% get wsindy_data object
traj_inds = [1:3];
% traj_inds = [4 5 6 7 8 11];
% traj_inds = [20 21]; %12-18, 22-26 too noisy, 19 too boxy
% traj_inds = [27 28 29 30 31 32]; %12-18 too noisy, 19 too boxy
% traj_inds = [20 2 5 27];

x_red = X_sigma_interp;
Uobj = arrayfun(@(i)wsindy_data(x_red{i},t_sigma_interp{i}),traj_inds(:));
arrayfun(@(U)U.coarsen(-500),Uobj);
nstates = Uobj.nstates;
M = Uobj.dims;
ntraj = length(Uobj);

figure(1);clf
for U=Uobj(:)'
    U.plotPhase;
    hold on
end
hold off

noise_ratio = 0;
rng('shuffle')
rng_seed = rng().Seed; rng(rng_seed);
Uobj.addnoise(noise_ratio,'seed',rng_seed);

%% get lib tags
polys = [0:5];
trigs = [];
tags_1 = get_tags(polys,[],nstates);
tags_2 = get_tags([],trigs,nstates);
tags = prodtags(tags_1,tags_2);
lib = library('tags',tags);

%% get test function
% tf = arrayfun(@(U)testfcn(U,'phifuns','delta','meth','direct','param',1,'mtmin',1,'subinds',2),Uobj,'uni',0);
tf = arrayfun(@(U)testfcn(U,'meth','FFT','param',1),Uobj,'uni',0);

%% build WSINDy linear system
WS = wsindy_model(Uobj,lib,tf);

%% solve
optm = WS_opt();
toggle_wendy = 3;
if toggle_wendy==0
    [WS,loss_wsindy,its,G,b] = optm.MSTLS(WS);
elseif toggle_wendy==1
    [WS,w_its,res,res_0,CovW] = optm.wendy(WS,'maxits',100,'regmeth','MSTLS');
elseif toggle_wendy==2
    [WS,loss_wsindy,its,G,b] = optm.MSTLS(WS,'alpha',0.01);
    [WS,w_its,res,res_0,CovW] = optm.wendy(WS,'maxits',20);
elseif toggle_wendy==3
    [WS,loss_wsindy,lambda,w_its,res,res_0,CovW] = optm.MSTLS_WENDy(WS,...
        'maxits_wendy',25,'lambda',10.^linspace(-4,-1,50),'verbose',1);
    disp(['wendy its at optimal lambda=',num2str(size(w_its,2))])
end

%% simulate learned and true reduced systems

toggle_compare = 1:ntraj;
if ~isempty(toggle_compare)
    w_plot = WS.weights;
    rhs_learned = WS.get_rhs;
    tol_dd = 10^-12;
    options_ode_sim = odeset('RelTol',tol_dd,'AbsTol',tol_dd*ones(1,nstates));

    for i=toggle_compare
        t_train = Uobj(i).grid{1};% linspace(tred(1),tred(end),2000);
        x0_reduced = Uobj(i).get_x0([]);
        [t_learned,xH0_learned{i}]=ode15s(@(t,x)rhs_learned(x),t_train,x0_reduced,options_ode_sim);
        figure(7+i-1);clf
        for j=1:nstates
            subplot(nstates,1,j)
            plot(Uobj(i).grid{1},Uobj(i).Uobs{j},'b-o',t_learned,xH0_learned{i}(:,j),'r-.','linewidth',2)
            try
                title(['rel err=',num2str(norm(xH0_learned{i}(:,j)-Uobj(i).Uobs{j})/norm(Uobj(i).Uobs{j}))])
            catch
            end
            legend({'data','learned'})
        end
    end
end

%%

figure(1); clf
semilogx(loss_wsindy(2,:),loss_wsindy(1,:),'o-')

%%

N = 50;
x = linspace(pi/3,7*pi/6,N);
y = linspace(-2*pi/3,pi/3,N);
[xx,yy] = meshgrid(x,y);

F = arrayfun(@(x,y)rhs_learned([x y]),xx,yy,'uni',0);
Fcell = arrayfun(@(j)cellfun(@(u)u(j),F),1:nstates,'un',0);
for j=1:nstates
    figure(j)
    colormap(bone)
    contourf(xx,yy,Fcell{j},50)
    colorbar
    hold on
    for i=1:ntraj
        plot(Uobj(i).Uobs{:})
        plot(Uobj(i).Uobs{:})
    end
    hold off
end