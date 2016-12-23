function [node_marg] = get_beliefs(G)
% GET_BELIEFS - Returns cell arrays containing beliefs for each node.
%
% Inputs:
%   G: Factor graph to perform loopy BP over
%
% Outputs:
%   node_marg: cell array containing the marginals of each variable node,
%              computed by multiplying and normalizing the current messages.
%              Same format as nodeMarg = marg_brute_force(G);
%
% Brown CS242

  num_var = numel(G.var);
  node_marg = cell(num_var,1);
  % FILL ME IN!
  for i=1:length(G.var)
      varmag = 1;
      for j=1:length(G.var(i).incoming)
          varmag = varmag .* G.var(i).incoming{j};
      end
      node_marg{i} = varmag/sum(varmag);      
  end
  
  
  

end
