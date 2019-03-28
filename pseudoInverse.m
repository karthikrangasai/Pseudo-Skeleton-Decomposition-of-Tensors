function [pseudoInv, U, S, V, r, UInverse, SInverse, VInverse, p, d] = pseudoInverse(A, tol)
%function [pseudoInv, U, S, V] = pseudoInverse(A, tol)
    
    s = max(size(A));

    [U S V] = svd(A);
    
    r = 0;
    
    for i=1:s
        if S(i,i) > tol
            %S(i,i) = 0;
            r = r + 1;
        end
    end
    
    UInverse = U(1:s,1:r);
    SInverse = S(1:r,1:r);
    VInverse = V(1:r,1:s);
    
    for i=1:r
        SInverse(i,i) = 1/S(i,i);
    end

    pseudoInv = VInverse * inv(SInverse) * UInverse';
    %pseudoInv = V * S * U';
    
    p = pinv(A, tol);
    d = max(max(abs(p - pseudoInv)));
end
