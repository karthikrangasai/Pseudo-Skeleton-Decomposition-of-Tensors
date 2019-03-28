function [d, A, k, c, C, r, R, pseudoInv, U, S, V, truncRank, UInverse, SInverse, VInverse, p, dPseudoInv, toler] = pseudoSkeletonDecomp(s, rank, tol)
%function [d, A, k, c, C, r, R, pseudoInv, U, S, V, toler] = pseudoSkeletonDecomp(s, rank, tol) %% Function call with lesser return values
%function [d, A, k, c, C, r, R, pseudoInv, toler] = pseudoSkeletonDecomp(s, rank, tol) %% %% Function call with more lesser return values

    x = 0:1/s:1;
    y = 0:1/s:1;
    A = zeros(s);
    
    for i=1:s
        for j=1:s
            A(i,j) = 1.0/(1.0 + x(i)^2 + y(j)^2);		%%%		Different functions f(x) 
            %A(i,j) = (x(i)^2 + y(j)^2);				%%%									for testing the					
            %A(i,j) = exp(x(i)^2 + y(j)^2);				%%%											    		algorithm on f(A)
        end
    end
    
    k = rank;
    c = sort(randperm(s,k));
    C = zeros(s,k);
    M = zeros(k,k);
    R = zeros(k,s);
    
    for i=1:k
        C(:,i) = A(:,c(i));
    end
    
    r = maxvol2(C);
    
    for i=1:k
        R(i,:) = A(r(i),:);
    end
    
    for i=1:k
        for j=1:k
            M(i,j) = A(r(i),c(j));
        end
    end
    
    % **** General Inverse ****
    %M_in = inv(M);
    
    % **** General Inverse ****   [ Testing between available pseudo inverse function and self written ]
    %Mps = pinv(M, 10^(-10));
    %Mps = pseudoInverse(M, 10^(-10));
    
    toler = tol; 
    
    [pseudoInv, U, S, V, r, UInverse, SInverse, VInverse, p, dPseudoInv] = pseudoInverse(M, tol);
    %[pseudoInv, U, S, V] = pseudoInverse(A, tol);
    
    %s = max(size(Mps));
    
    truncRank = r;
    
    %T = C * Mps * R;
    T = C * pseudoInv * R;
    
    %d = norm((A - T), inf);
    d = max(max(abs(A - T)));
end
