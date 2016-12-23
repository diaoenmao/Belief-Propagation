function nodeMarg = marg_brute_force( G )
% MARGE_BRUTE_FORCE - Compute marginals by brute force enumeration.  This
%   is mainly used for debugging on smaller graphical models.
%
% Brown CS242

  num_vars = numel(G.var);

  % compute full joint
  dims = [ G.var.dim ];
  logP = zeros(dims);
  for fac_i = 1:numel(G.fac)
    fac = G.fac(fac_i);
    nbrs = fac.nbrs_var;
    
    % permute dimensions
    I = [ nbrs, setdiff(1:num_vars, nbrs) ];
    vardims_perm = dims(I);
    logP_perm = permute(logP, I);
    
    % multiply factor into joint
    vardims_perm(1:numel(nbrs)) = 1;
    fac_rep = repmat(fac.p, vardims_perm);
    logP_perm = logP_perm + log(fac_rep);
    
    % unpermute
    logP = ipermute(logP_perm, I);    
  end
  logP = logP - max(logP(:));
  P = exp(logP) ./ sum(exp(logP(:)));
  
  % compute marginals
  for dim=1:num_vars
    % special case: observed variables
    if G.var(dim).observed
      thisMarg = zeros(dims(dim),1);
      thisMarg(G.var(dim).observed) = 1.0;
    else
      marg_vars = 1:num_vars;
      marg_vars(dim) = [];
      marg_vars = sort(marg_vars, 'descend');
      
      % marginalize
      thisMarg = P;
      for dim_marg = marg_vars
        thisMarg = sum(thisMarg, dim_marg);
      end
      thisMarg = reshape(thisMarg, [dims(dim), 1]);
    end

    nodeMarg{dim} = thisMarg;
  end
end

