function [x, res1, res2] = funIFVAMP(A, y, img, D, damping, lam, Rho, Lip, Niter, N, r)
    n = sqrt(N);
    x = zeros(N,Niter);
    u = zeros(r,Niter);
    z = zeros(r,Niter);
    v = zeros(N,Niter);
    t = zeros(N,Niter);
    rho = zeros(1,Niter);  % timestep
    sgmx = zeros(1,Niter); % variance of x
    sgmz = zeros(1,Niter); % variance of z
    res1 = zeros(Niter,1);  % save SSIM per iteration
    res2 = zeros(Niter,1);  % save RHO updated per iteration
    L1 = Lip; 
    L2 = Lip;
    rho(1) = Rho; im=img; lambda = lam; res2(1) = Rho;
    for i = 1:Niter
        % update x
        tempinv = diag(1./diag((1/2*L2 + 1/2*rho(i)*L1)*eye(N))) ;
        x(:,i) = tempinv * (   A'*(y - A*t(:,i)) + ...
                               1/2* L2 *t(:,i) + ...                        
                               funKT( (u(:,i) - rho(i)*funK(v(:,i))),D ) + ...
                               1/2*rho(i)* L1 *v(:,i) ); 
        sgmx(i) = trace(tempinv)./N;
        % update z 
        p1 = ( funK(x(:,i))-sgmx(i)*u(:,i) )./( 1-sgmx(i)*rho(i) );
        p2 = ( lambda*sgmx(i) )./( 1-sgmx(i)*rho(i) );
        z(:,i) = ( 1>p2./norm(p1) )*( 1-p2./norm(p1) ).* p1; 
        sgmz(i) = (sgmx(i))./(1-sgmx(i)*rho(i)) * (norm(p1)>p2)*(r-(r-1)*p2./norm(p1))./r; 
        u(:,i+1) = u(:,i) + damping * ( z(:,i)./sgmz(i) - funK(x(:,i))./sgmx(i) );
        rho(i+1) = rho(i) + damping * (1./sgmz(i) - 1./sgmx(i));
        res1(i) = ssim(abs(255.*reshape(x(:,i),[n n])),255.*im);
        res2(i+1) = rho(i+1);
        % update L 
        L1 = L1track( u(:,i+1),rho(i+1),v(:,i),x(:,i),D );
        L2 = L2track( y,A,x(:,i),t(:,i) );
        % update v
        v(:,i+1) = x(:,i);
        t(:,i+1) = x(:,i);
    end
end

function minLv = L1track(u,ro,v,x,D)
    fval = norm(funK(x)-u./ro)^2;
    gpart = norm(funK(v)-u./ro)^2 + 2*(x-v)'*funKT(funK(v)-u./ro,D);
    minLv = 2*(fval - gpart)./norm(x-v)^2 + 0.5;
end

function L2 = L2track(y,A,x,t)
    fval = norm(A*x-y)^2;
    gpart = norm(A*t-y)^2 + 2*(x-t)'*A'*(A*t-y);
    L2 = 2*(fval - gpart)./norm(x-t)^2 + 0.5;
end



