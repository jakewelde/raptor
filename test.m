time = zeros(100,1);
for i=1:size(time,1)
    tic
    compute_ff;
    time(i) = toc;
end
mean(time)
variance(time)