clear
tol = 1e-8; maxiter = 5000;

% Imprime en archivo para tener formato de tabla de Latex
fidw_latex= fopen('Resultados_Poisson_Latex.txt','w');
fprintf(fidw_latex, 'Matriz & n & iter & $\\norm{x_{GC} - x_{GCP}}_2$ \\\\ \\hline \n');

% Imprime en consola
fprintf('%12s %6s %9s %12s\n','Matriz','n','iter', '||x_GC - x_GCPre||');

for q = [10, 20, 30];

    n = q*q;

    A = gallery('poisson',q);
    x0 = zeros(n,1);
    b = A*ones(n,1);

    %Prueba del GC
    M = eye(n);
    [x_GC, iter_GC, t_GC] = GCPre(x0, A, b, tol, maxiter, M);

    %Prueba del GC Precondicionado con Cholesky incompleto
    opts.michol = 'on';
    C = ichol(A);
    M = C'*C;
    [x_GCP_ichol, iter_GCP_ichol, t_GCP_ichol] = GCPre(x0, A, b, tol, maxiter, M);

    %Prueba del GC Precondicionado con Jacobi
    M = diag(diag(A));
    [x_GCP_jac, iter_GCP_jac, t_GCP_jac] = GCPre(x0, A, b, tol, maxiter, M);
    
    % Norma de las diferencias
    norma_dif_ichol = norm(x_GC - x_GCP_ichol);
    norma_dif_jac = norm(x_GC - x_GCP_jac);

    fprintf('%12s %6d %9d %12e \n','Cholesky', n, iter_GCP_ichol, norma_dif_ichol);
    fprintf('%12s %6d %9d %12e \n','Jacobi', n, iter_GCP_jac, norma_dif_jac);

    % Formato tabla de Latex

    fprintf(fidw_latex, '%s & %d & %d & %e \\\\ \n', 'Cholesky', n, iter_GCP_ichol, norma_dif_ichol);
    fprintf(fidw_latex, '%s & %d & %d & %e \\\\ \n', 'Jacobi', n, iter_GCP_jac, norma_dif_jac);

end