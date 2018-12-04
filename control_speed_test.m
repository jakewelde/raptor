time = linspace(0,total_dt,100);
results = zeros(size(time));
for i=1:size(time,2)
    tic;
	compute_control(stacked,segment_dt*(i-1),1);
    results(i) = toc;
end
mean(results)