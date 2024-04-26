%%% get true dynamics
function [w_true,rhs_H_true,Hfun_true] = get_true_dyn_script(tags,true_prod_tags,true_coeff,params,rhs_p)
    J = length(tags);
    w_true = zeros(J,1);
    try
        true_inds = gettrueinds(tags,true_prod_tags);
        w_true(true_inds)=true_coeff(:);
    end
    rhs_H_true = @(varargin) varargin{1}*0;
    try
        rhs_H_true = @(x)rhs_p(x,params);
    end
    Hfun_true = @(varargin) varargin{1}*0;
    try
        libtrue = library('tags',true_prod_tags);
        Hlib_true = cellfun(@(tm) tm.fHandle, libtrue.terms, 'uni',0);
        Hfun_true = @(varargin) varargin{1}*0;
        for i=1:length(true_coeff)
            Hfun_true = @(varargin) Hfun_true(varargin{:}) + true_coeff(i)*Hlib_true{i}(varargin{:});
        end    
    end
end
