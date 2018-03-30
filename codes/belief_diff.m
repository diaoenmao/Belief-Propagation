function D = belief_diff( b1, b2 )
% BELIEF_DIFF - Computes the symmetric L1 distance between belief's 'b1'
%   and 'b2'. Inputs are N-dim cell arrays, where each cell is a probability
%   vector.  L1 distances are returned in an Nx1 vector.
%
% Brown CS242

  num_b = numel(b1);
  if num_b ~= numel(b2), error('belief_diff: Belief dimensions are not equal.\n'); end;
  
  D = zeros(num_b,1);
  for i=1:num_b
    D(i) = 1/length(b1{i}) * sum( abs( b1{i} - b2{i} ) );
  end
  D = mean(D);
end

