function outputMatrix = nta_neighsearch(inputMatrix, subq, mode)
    [n, ~] = size(inputMatrix);
    n_sub = length(subq);
    
    outputMatrix = zeros(n, n_sub);
    
    if mode == 1
        outputMatrix = pdist2(inputMatrix, inputMatrix(subq, :), 'euclidean');
        
    elseif mode == 2
        outputMatrix = pdist2(inputMatrix, inputMatrix(subq, :), 'chebychev');
    end
end