%% construction code test
xfun = @(t)[sin(t) cos(t)];
xdotfun = @(t)[cos(t) -sin(t)];
t = linspace(0,10,1000)';
x = cell2mat(arrayfun(@(t)xfun(t),t,'un',0));
xdot = cell2mat(arrayfun(@(t)xdotfun(t),t,'un',0));

Jfuns = {@(x)x(1).^2,@(x)x(2).^2,@(x)x(1).*x(2)};
Jgradfuns = {@(x)[2*x(1) 0],@(x)[0 2*x(2)],@(x)[x(2) x(1)]};
dfuns = {@(x)1+0.2*sin(x(1)),@(x)1-0.2*cos(x(2))};
Afuns = {@(x)1+0.2*sin(x(1)),@(x)1-0.2*cos(x(2));@(x)1+0.2*sin(x(1)),@(x)1-0.2*cos(x(2))};
G = buildG(x,xdot,Jgradfuns,dfuns,Afuns);

[U,S,V] = svd(G,0);

%% ode test
tol = 10^-12;
M = 1000;
t = linspace(0,300,M);
x0 = [1 1];
w = 2;
a = 0.2; c = 0.15;
options_ode_sim = odeset('RelTol',tol,'AbsTol',tol*ones(length(x0),1));
[t,x]=ode15s(@(t,x)rhs(x,w,a,c),t,x0,options_ode_sim);

xdot = (diff(x(1:end-1,:))+diff(x(2:end,:)))/2/mean(diff(t));
x = x(2:end-1,:);
t = t(2:end-1);

M = length(t);

Jfuns = {@(x)x(1).^2,@(x)x(2).^2,@(x)x(1).*x(2)};
Jgradfuns = {@(x)[2*x(1) 0],@(x)[0 2*x(2)],@(x)[x(2) x(1)]};
dfuns = {@(x)x(1),@(x)x(2)};
Afuns = {@(x)0*x(1),@(x)a+c*sin(w*x(1));@(x)-a-c*sin(w*x(1)),@(x)0*x(1)};
[G,Jcell,denom,Acell,P] = buildG(x,xdot,Jgradfuns,dfuns,Afuns);
Jmat = cell2mat(cellfun(@(f)arrayfun(@(i)f(x(i,:)),(1:M)'),Jfuns,'un',0));

subsamp = 0.2;
[U,S,V] = svd(G(randperm(size(G,1),ceil(size(G,1)*0.1)),:),0);
W = V(:,end);
disp(W')

rhs_learned = @(x) rhs_J(x,Jgradfuns,W,dfuns,Afuns,0);
options_ode_sim = odeset('RelTol',tol,'AbsTol',tol*ones(length(x0),1));
[t_L,x_L]=ode15s(@(t,x)rhs_learned(x),t,x0,options_ode_sim);

clf; plot(t,x,t_L,x_L,'--','linewidth',3)

%% functions 

function dx = rhs(x,w,a,c)
%%% J=1/2(x^2+y^2), A = [0 sin(w*x);-sin(w*x) 0], d = [x y];
   d = x(1)^2+x(2)^2;
   dx = x*0;
   dx(1) = (a+c*sin(w*x(1)))*x(2)/d;
   dx(2) = (-a-c*sin(w*x(1)))*x(1)/d;
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
