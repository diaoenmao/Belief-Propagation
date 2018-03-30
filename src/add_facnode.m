function G = add_facnode(G, p, varargin)
% ADD_FACNODE - Add factor node to factor graph 'G'.  
%
% INPUTS:
%   G - Factor graph
%
%   p - Potential matrix with p(i,j,...) the potential for
%       x_a=i, x_b=j, ...
%
%   varargin - Variable nodes involved in factor.  Order matches dimensions
%              potential, e.g. varargin = { a, b, ... } => p(x_a,x_b,...).
%
% Brown CS242

  f.p = p;
  f.nbrs_var = [ varargin{:} ];
  f.id = numel(G.fac) + 1;
  
  % check dimensions
  if ( ( numel(f.nbrs_var) > 1 ) && ( numel(f.nbrs_var) ~= ndims(p)) ) || ...
    ( ( numel(f.nbrs_var) == 1 ) && ( ndims(p) ~= 2 ) )
    error('add_facnode: Factor dimensions does not match size of domain.');
  end
  
  G.fac = cat( 1, G.fac, f );
  for I=f.nbrs_var
    G.var(I).nbrs_fac = [ G.var(I).nbrs_fac; f.id ];
  end
end
