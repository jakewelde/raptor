function handle_array = sample(constraint,t_range)
    handle_array = cell(size(t_range));
    for i = 1:size(t_range,2)
        handle_array{i} = @(C) constraint(C,t_range(i));
    end
end