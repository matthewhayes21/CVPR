function [eigenVectors,eigenValues] = LDA(Data1, Data2)
    % returnes the eigenVectors and the eigenValues for LDA for data with
    % two classes
    data = [transpose(makeVector(Data1)); transpose(makeVector(Data2))];
    As = [mean(Data1(:, 1)),mean(Data1(:, 2));
          mean(Data2(:, 1)),mean(Data2(:, 2))];

    % calculate Sw
    D = data - [ kron(ones(1,10), As(:, 1)) , kron(ones(1,10), As(:, 2)) ];
    Sw = D*transpose(D);

    % calculate Sb
    meanDiff = (As(:, 1) - As(:, 2));
    Sb = meanDiff * transpose(meanDiff);

    [eigenVectors, eigenValues] = eig(Sw\Sb);
end