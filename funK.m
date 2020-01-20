function Kx = funK(x)
% function out=tlv(X,type)
    n = sqrt(length(x));
    x = reshape(x,[n n]);
    P1 = x(1:n-1,:)-x(2:n,:);
    P2 = x(:,1:n-1)-x(:,2:n);
    Kx = -[P1(:); P2(:)];
end








