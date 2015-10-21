function [x, iter, t] = GCPre(x, A, b, tol, maxiter, M)
% Resuelve el sistema lineal Ax = b utilizando gradiente conjugado
% precondicionado, utiliza el vector x como semilla
% In
% x.- Vector inicial
% A.- Matriz s.p.d. n*n
% b.- vector de tamaño n
% tol.- tolerancia para la norma del residual
% maxiter.- Número máximo de iteraciones
% M.- Precondicionador
% Out
% x.- aproximacion a la solución del sistema Ax=b
% iter.- Número total de iteraciones
    
    tic
    % x = zeros(length(A), 1);
    r0 = A * x - b;
    u = norm(r0);
    v = Inf;
    y0 = M\r0;
    %y0 = M'\y0;
    p = -y0;
    iter = 0;
    while v > tol * u && iter < maxiter
        a = (r0' * y0) / (p' * A * p);
        x = x + a * p;
        r = r0 + a * A * p;
        y = M\r;
        y = M'\y;
        beta = (r' * y)/(r0' * y0);
        p = - y + beta * p;
        r0 = r;
        y0 = y;
        v = norm(r);
        iter = iter + 1;
    end
    t = toc;
end