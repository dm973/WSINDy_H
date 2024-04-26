function rhs_H = get_rhs_H(Hlib,w,nstates,J)
    if ~exist('J','var')
        J = [];
    end
    if any(w)
        grads_H = grad_sym(Hlib(w~=0),nstates);
        rhs_H = @(x)rhs_dH(grads_H,w(w~=0),x,J);
    else
        rhs_H = @(x)x*0;
    end
end

function dH = rhs_dH(grads,W,x,J)
    if isempty(J)
        J = kron(eye(nstates/2),[[0 1];[-1 0]]);
    end
    nstates = length(x);
    x = num2cell(x);
    gradH = sum(cell2mat(cellfun(@(f,w)w*f(x{:}),grads(:)',num2cell(W(:)'),'uni',0)),2);
    dH = J*gradH;
end