t_f = 1;
C = find_coefficients([0;0;0;0],[0;-1;0;0],t_f);

derivs = 4;
time_range = linspace(0,t_f,100);
derivatives = zeros(derivs,size(time_range,2));
for i = 1:size(time_range,2)
    t_dim = time_range(i);
    d = compute_derivatives(C,t_dim,t_f);
    derivatives(:,i) = d(1:derivs);
end

plot(time_range,derivatives);
shg;
legend('position','velocity','acceleration','jerk');
[derivatives(:,1) derivatives(:,end)]