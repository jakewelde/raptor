syms time time_f
basis = time.^(0:7);
derivative_matrix = [
    basis;
    diff(basis,time);
    diff(basis,time,2);
    diff(basis,time,3);
    diff(basis,time,4);
];
D = [
    subs(derivative_matrix(1:4,:),time,0);
    subs(derivative_matrix(1:4,:),time,1);
];
D = double(inv(D));

derivative_function = matlabFunction(derivative_matrix);

dimensionalize = @(t_f,vector) diag((1/t_f).^(0:size(vector,1)-1))*vector;
nondimensionalize = @(t_f,vector) diag((t_f).^(0:size(vector,1)-1))*vector;

global compute_derivatives
compute_derivatives = @(c,t,t_f) dimensionalize(t_f,derivative_function(t/t_f)*c);

find_coefficients = @(init_state,final_state,t_f) D*[nondimensionalize(t_f,init_state); nondimensionalize(t_f,final_state)];
% find_coefficients_intermediate = @(init_state,crit_state,t_crit,t_f) double([
%     subs(derivative_matrix(1:4,:),time,0);
%     subs(derivative_matrix(1:4,:),time,t_crit);
% ]*[nondimensionalize(t_f,init_state); nondimensionalize(t_f,crit_state)]);