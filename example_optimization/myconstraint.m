function [val] = myconstraint(x,a,b)
    val = (a*x(1) + b)^3 - x(2);
end