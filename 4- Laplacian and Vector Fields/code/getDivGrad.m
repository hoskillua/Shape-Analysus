%% PROBLEM 2 - YOUR CODE HERE TO COMPUTE DIVERGENCE AND GRADIENT
% Div (|V| x 3|F| sparse) maps face-based vector fields to their vertex-based divergence
% Grad (3|F| x |V| sparse) maps a scalar function to its face-based gradient
function [Div, Grad] = getDivGrad(data)

Div = sparse(data.nv, 3 * data.nf);
Grad = sparse(3 * data.nf, data.nv);

end
%%% END HOMEWORK PROBLEM