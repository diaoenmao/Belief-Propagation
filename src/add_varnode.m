function [G, id] = add_varnode( G, name, dim )
% ADD_VARNODE - Add variable node to factor graph 'G'.
%
% Brown CS242

  id = numel(G.var) + 1;
  
  v.name = name;
%   v.dim = numel(vals);
%   v.vals = vals;
  v.dim = dim;
  v.id = id;
  v.nbrs_fac = [];
  v.observed = 0;
  G.var = cat( 1, G.var, v );
end
