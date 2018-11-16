function u = op(matrix_column,w)
%     [U, ~, V] = svd(reshape(matrix_column,[3 3]));
%     u = (U * V') * w;
    u = reshape(matrix_column,[3 3]) * w;
end

