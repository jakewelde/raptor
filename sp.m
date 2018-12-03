function u = sp(matrix_column,w)
    R = reshape(matrix_column,[3 3]);
    [U,~,V] = svd(R);
    R = U*V';
    u = R * w;
end

