syms time
basis = time.^(0:7);
derivative_matrix = [
    basis;
    diff(basis,time)
    diff(basis,time,2)
    diff(basis,time,3)
    diff(basis,time,4)
];
D = [
    subs(derivative_matrix(1:4,:),time,0);
    subs(derivative_matrix(1:4,:),time,1);
];
D = double(inv(D));

derivative_matrix = matlabFunction(derivative_matrix);

compute_derivatives = @(c,t) derivative_matrix(t)*c;