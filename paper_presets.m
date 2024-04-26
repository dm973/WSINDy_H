if expnum==1
    mQ = 1;
    ep = 0.05;
    train_inds = 14;
    polys = 1:2;
    trigs = 0:2:6; L = 1;
    add_polytrig = 1;
    pdeg_x0 = 0;
elseif expnum==2
    train_inds = 27;
    mQ = 5;
    ep = 0.05;
    polys = 1:4;
    trigs = []; L = 1;
    add_polytrig = 0;
    pdeg_x0 = 0;
elseif expnum==3
    mQ=1;
    ep = 0.03;
    train_inds = 9;%7;
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
elseif expnum==5
    mQ = 1;
    mu_ep = 0.5; % adiabatic invariant, expnum 5 only
    ep = 0.03;
    train_inds = 3;
    polys = 1:2;
    trigs = 0:3; L = 1;
    add_polytrig = 0;
    pdeg_x0 = 3;    
end