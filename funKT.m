function KTz = funKT(z,D)
% transposition of the two dimensional TV matrix
    % to decompose the large scale KT matrix-vector mutiplication:
    %   sub-block1: 1st
    %   sub-block2: from 2nd to (n-1)-th
    %   sub-block3: n-th
    n = size(D,2); 
    N = n*n;
    r = 2*(n-1)*n; % length of vector z
    block1 = zeros(n,r);
    block3 = zeros(n,r);
    % construction of each block
    block1(1:n,1:(n-1)) = D';
    block1(1:n,N-(n-1):N) = -1.*eye(n);
    KTz(1:n,1) = block1 * z; 
    for c=2:n-1
        block2 = zeros(n,r);
        gap = (N+c-1)-(n-1)-2*n;      
        b1s = (c-1)*(n-1)+1;
        b1e = b1s+(n-2);
        b2s = b1e+1+gap;
        b2e = b2s+(n-1);
        b3s = b2e+1;
        b3e = b3s+(n-1);
        block2(1:n, b1s:b1e) = D';
        block2(1:n, b2s:b2e) = eye(n);
        block2(1:n, b3s:b3e) = -1.* eye(n);
        lis = (c-1)*n+1;
        lie = lis+n-1;
        KTz(lis:lie) = block2 * z; 
    end
    b1s = c*(n-1)+1;
    b1e = b1s+(n-2);
    b2s = b1e+1+gap+1;
    b2e = b2s+(n-1);
    block3(1:n, b1s:b1e) = D';
    block3(1:n, b2s:b2e) = eye(n);
    KTz(N-(n-1):N,1) = block3 * z; % last n entries 
end








