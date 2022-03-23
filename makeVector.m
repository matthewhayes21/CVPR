function vec = makeVector(data)
    % takes in a mxN matrix and stacks all columns to produce a vector
    vec = [];
    N = width(data);
    for i = 1:N
        vec = [vec; data(:, i)];
    end
end