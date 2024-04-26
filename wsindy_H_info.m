%% load data
clc;
Info = zeros(24,1);
kappa = Info;
for kk=1:24

expnum = 4;
mQ = 2;
ep = 0.01;
train_inds = kk;
test_inds = [train_inds];
test_ratio = 1;
toggle_plot_sol = 0;
ppTf = -8; % points per fast mode
ncyc = -1; % number of slow modes observed
noise_ratio = 0; % noise ratio

load_data_script;

%%% get lib

polys = 1:2;
trigs = 1:3;
include_cos = 1;

get_lib_script;

%%% get initial test functions

phifun = optTFcos(2,2);
mtmax_fac = 0.5;
sobol_inds = 0;
tf_meth = 'direct';
tf_param = floor(Uobj{1}.dims/20);

tf_script;

%%% solve linear system

WS = wsindy_model(repmat(Uobj,length(sobol_inds),1),lib,tf,'toggleH',1);
WS.cat_Gb;
b = WS.b{1};
[w_true,rhs_H_true,Hfun_true] = get_true_dyn_script(tags,true_prod_tags,true_coeff,params,rhs_p);

s = find(w_true);
Itemp = zeros(length(s),1);
kappa(kk) = cond(WS.G{1});
for j=1:length(s)
    G = WS.G{1}(:,[1:s(j)-1 s(j)+1:end]);
    Itemp(j) = norm(b-G*(G\b))/norm(b);
end
Info(kk) = min(Itemp);

end
