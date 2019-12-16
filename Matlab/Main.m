[first_var, second_var] = first_function(4, 6)
load (variables.mat, variables)

function [a, b] = first_function(x, y) 
    a = x + y;
    b = x * y;
end
