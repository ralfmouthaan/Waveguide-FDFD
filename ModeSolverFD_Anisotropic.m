% Ralf Mouthaan
% University of Cambridge
% May 2020
%
% FD Mode solver
%
% As per:
% Zhu & Brown, "Full-Vectorial Finite-Difference Analysis of 
% Microstructured Optical Fibers"
% and
% Yu & Chang, "Yee-Mesh-Based Finite Difference Eigenmode Solver with PML
% Absorbing Boundary Conditions for Optical Waveguides and Photonic 
% Crystal Fibers
%
% Notes:
% * Uses a Yee mesh
% * Zero-value boundary conditions
% * Stretched-Coordinate PML
% * Everything is stored in row-major format

function [RetVal] = ModeSolverFD_Anisotropic(dx, nx, ny, nz, lambda, beta, NoModes)

    %% Error-checking

    if size(nx, 1) ~= size(nx, 2)
        error('Expecting square problem space');
    end
    if size(nx) ~= size(ny)
        error('Expecting size(nx) = size(ny)');
    end
    if size(ny) ~= size(nz)
        error('Expecting size(ny) = size(nz)');
    end
    if lambda/dx < 10
        warning('lambda/dx < 10: this will likely cause discretisation errors');
    end

    %% User-Defined

    eps0 = 8.85e-12;
    mu0 = 4*pi*10^(-7);
    c = 3e8;
    
    if ~exist('NoModes', 'var')
        NoModes = 1;
    end

    Nx = size(nx, 1);
    f = c/lambda;
    w = 2*pi*f;
    k0 = 2*pi/lambda;

    PML_Depth = 30;
    PML_TargetLoss = 1e-8;
    PML_PolyDegree = 2;
    
    PML_SigmaMax = (PML_PolyDegree + 1)/2 * eps0*c/PML_Depth/dx* log(1/PML_TargetLoss); % Ballpark, have not accounted for refractive index.

    Epsrx = nx.^2;
    Epsry = ny.^2;
    Epsrz = nz.^2;
    Epsrx = MatrixToColumn(Epsrx);
    Epsry = MatrixToColumn(Epsry);
    Epsrz = MatrixToColumn(Epsrz);

    %% Ux, Uy, Vx, Vy
    
    fprintf('Calculating Ux, Uy, Vx, Vy...\n');

    I = speye(Nx*Nx, Nx*Nx);

    idx_x = 0;
    idx_y = 1;

    Epsx = zeros(1,Nx*Nx);
    Epsy = zeros(1,Nx*Nx);
    Epsz = zeros(1,Nx*Nx);

    Ax_idxi = zeros(1,2*Nx*Nx);
    Ax_idxj = zeros(1,2*Nx*Nx);
    Ax_vals = zeros(1,2*Nx*Nx);

    Ay_idxi = zeros(1,2*Nx*Nx);
    Ay_idxj = zeros(1,2*Nx*Nx);
    Ay_vals = zeros(1,2*Nx*Nx);

    Bx_idxi = zeros(1,2*Nx*Nx);
    Bx_idxj = zeros(1,2*Nx*Nx);
    Bx_vals = zeros(1,2*Nx*Nx);

    By_idxi = zeros(1,2*Nx*Nx);
    By_idxj = zeros(1,2*Nx*Nx);
    By_vals = zeros(1,2*Nx*Nx);

    Cx_idxi = zeros(1,2*Nx*Nx);
    Cx_idxj = zeros(1,2*Nx*Nx);
    Cx_vals = zeros(1,2*Nx*Nx);

    Cy_idxi = zeros(1,2*Nx*Nx);
    Cy_idxj = zeros(1,2*Nx*Nx);
    Cy_vals = zeros(1,2*Nx*Nx);

    Dx_idxi = zeros(1,2*Nx*Nx);
    Dx_idxj = zeros(1,2*Nx*Nx);
    Dx_vals = zeros(1,2*Nx*Nx);

    Dy_idxi = zeros(1,2*Nx*Nx);
    Dy_idxj = zeros(1,2*Nx*Nx);
    Dy_vals = zeros(1,2*Nx*Nx);

    for i = 1:Nx*Nx

        % Keep track of where we are
        idx_x = idx_x + 1;
        if idx_x > Nx
            idx_y = idx_y + 1;
            idx_x = 1;
        end
        West_Dist = idx_x-1;
        North_Dist = idx_y-1;
        East_Dist = Nx - idx_x;
        South_Dist = Nx - idx_y;

        % Epsx, Epsy, Epsz
        if i - Nx >1
            Epsx(i) = (Epsrx(i) + Epsrx(i-Nx))/2;
        else
            Epsx(i) = Epsrx(i);
        end
        if i - 1 > 1
            Epsy(i) = (Epsry(i) + Epsry(i-1))/2;
        else
            Epsy(i) = Epsry(i);
        end
        if i - 1 - Nx > 1
            Epsz(i) = (Epsrz(i) + Epsrz(i-1) + Epsrz(i-Nx) + Epsrz(i-1-Nx))/4;
        else
            Epsz(i) = Epsrz(i);
        end

        % Sx, Sy
        if West_Dist <= PML_Depth
            Sx_Ey = 1-PML_SigmaMax*(1-West_Dist/PML_Depth)^PML_PolyDegree*1i/w/eps0/sqrt(Epsy(i));
            Sx_Ez = 1-PML_SigmaMax*(1-West_Dist/PML_Depth)^PML_PolyDegree*1i/w/eps0/sqrt(Epsz(i));
            Sx_Hy = 1-PML_SigmaMax*(1-(West_Dist-0.5)/PML_Depth)^PML_PolyDegree*1i/w/eps0/sqrt(Epsx(i));
            Sx_Hz = 1-PML_SigmaMax*(1-(West_Dist-0.5)/PML_Depth)^PML_PolyDegree*1i/w/eps0/sqrt(Epsrz(i)); % Is this Epsrz correct?
        elseif East_Dist <= PML_Depth
            Sx_Ey = 1-PML_SigmaMax*(1-(East_Dist-0.5)/PML_Depth)^PML_PolyDegree*1i/w/eps0/sqrt(Epsy(i));
            Sx_Ez = 1-PML_SigmaMax*(1-(East_Dist-0.5)/PML_Depth)^PML_PolyDegree*1i/w/eps0/sqrt(Epsz(i));
            Sx_Hy = 1-PML_SigmaMax*(1-East_Dist/PML_Depth)^PML_PolyDegree*1i/w/eps0/sqrt(Epsx(i));
            Sx_Hz = 1-PML_SigmaMax*(1-East_Dist/PML_Depth)^PML_PolyDegree*1i/w/eps0/sqrt(Epsrz(i)); % Is this Epsrz correct?
        end
        if North_Dist <= PML_Depth
            Sy_Ex = 1-PML_SigmaMax*(1-North_Dist/PML_Depth)^PML_PolyDegree*1i/w/eps0/sqrt(Epsx(i));
            Sy_Ez = 1-PML_SigmaMax*(1-North_Dist/PML_Depth)^PML_PolyDegree*1i/w/eps0/sqrt(Epsz(i));
            Sy_Hx = 1-PML_SigmaMax*(1-(North_Dist-0.5)/PML_Depth)^PML_PolyDegree*1i/w/eps0/sqrt(Epsy(i));
            Sy_Hz = 1-PML_SigmaMax*(1-(North_Dist-0.5)/PML_Depth)^PML_PolyDegree*1i/w/eps0/sqrt(Epsrz(i)); % Is this Epsrz correct?
        elseif South_Dist <= PML_Depth
            Sy_Ex = 1-PML_SigmaMax*(1-(South_Dist-0.5)/PML_Depth)^PML_PolyDegree*1i/w/eps0/sqrt(Epsx(i));
            Sy_Ez = 1-PML_SigmaMax*(1-(South_Dist-0.5)/PML_Depth)^PML_PolyDegree*1i/w/eps0/sqrt(Epsz(i));
            Sy_Hx = 1-PML_SigmaMax*(1-South_Dist/PML_Depth)^PML_PolyDegree*1i/w/eps0/sqrt(Epsy(i));
            Sy_Hz = 1-PML_SigmaMax*(1-South_Dist/PML_Depth)^PML_PolyDegree*1i/w/eps0/sqrt(Epsrz(i)); % Is this Epsrz correct?
        end

        %Ax
        Ax_idxi(2*i) = i;
        Ax_idxj(2*i) = i;
        Ax_vals(2*i) = -1/Sx_Ez;
        if i+1 <= Nx*Nx
            Ax_idxi(2*i+1) = i;
            Ax_idxj(2*i+1) = i+1;
            Ax_vals(2*i+1) = 1/Sx_Ez;
        end

        %Bx
        Bx_idxi(2*i) = i;
        Bx_idxj(2*i) = i;
        Bx_vals(2*i) = -1/Sx_Ey;
        if i+1 <= Nx*Nx
            Bx_idxi(2*i+1) = i;
            Bx_idxj(2*i+1) = i+1;
            Bx_vals(2*i+1) = 1/Sx_Ey;
        end

        %Ay
        Ay_idxi(2*i) = i;
        Ay_idxj(2*i) = i;
        Ay_vals(2*i) = -1/Sy_Ez;
        if i+Nx <= Nx*Nx
            Ay_idxi(2*i+1) = i;
            Ay_idxj(2*i+1) = i+Nx;
            Ay_vals(2*i+1) = 1/Sy_Ez;
        end

        %By
        By_idxi(2*i) = i;
        By_idxj(2*i) = i;
        By_vals(2*i) = -1/Sy_Ex;
        if i+Nx <= Nx*Nx
            By_idxi(2*i+1) = i;
            By_idxj(2*i+1) = i + Nx;
            By_vals(2*i+1) = 1/Sy_Ex;
        end

        %Cx
        Cx_idxi(2*i) = i;
        Cx_idxj(2*i) = i;
        Cx_vals(2*i) = 1/Sx_Hz;
        if i-1 >= 1
            Cx_idxi(2*i+1) = i;
            Cx_idxj(2*i+1) = i-1;
            Cx_vals(2*i+1) = -1/Sx_Hz;
        end

        %Dx
        Dx_idxi(2*i) = i;
        Dx_idxj(2*i) = i;
        Dx_vals(2*i) = 1/Sx_Hy;
        if i-1 >= 1
            Dx_idxi(2*i+1) = i;
            Dx_idxj(2*i+1) = i-1;
            Dx_vals(2*i+1) = -1/Sx_Hy;
        end

        %Cy
        Cy_idxi(2*i) = i;
        Cy_idxj(2*i) = i;
        Cy_vals(2*i) = 1/Sy_Hz;
        if i-Nx >= 1
            Cy_idxi(2*i+1) = i;
            Cy_idxj(2*i+1) = i-Nx;
            Cy_vals(2*i+1) = -1/Sy_Hz;
        end

        %Dy
        Dy_idxi(2*i) = i;
        Dy_idxj(2*i) = i;
        Dy_vals(2*i) = 1/Sy_Hx;
        if i-Nx >= 1
            Dy_idxi(2*i+1) = i;
            Dy_idxj(2*i+1) = i-Nx;
            Dy_vals(2*i+1) = -1/Sy_Hx;
        end

    end

    Ax_vals(Ax_idxi == 0) = [];
    Ax_idxj(Ax_idxi == 0) = [];
    Ax_idxi(Ax_idxi == 0) = [];
    Ax = sparse(Ax_idxi, Ax_idxj, Ax_vals, Nx*Nx, Nx*Nx);

    Ay_vals(Ay_idxi == 0) = [];
    Ay_idxj(Ay_idxi == 0) = [];
    Ay_idxi(Ay_idxi == 0) = [];
    Ay = sparse(Ay_idxi, Ay_idxj, Ay_vals, Nx*Nx, Nx*Nx);

    Bx_vals(Bx_idxi == 0) = [];
    Bx_idxj(Bx_idxi == 0) = [];
    Bx_idxi(Bx_idxi == 0) = [];
    Bx = sparse(Bx_idxi, Bx_idxj, Bx_vals, Nx*Nx, Nx*Nx);

    By_vals(By_idxi == 0) = [];
    By_idxj(By_idxi == 0) = [];
    By_idxi(By_idxi == 0) = [];
    By = sparse(By_idxi, By_idxj, By_vals, Nx*Nx, Nx*Nx);

    Cx_vals(Cx_idxi == 0) = [];
    Cx_idxj(Cx_idxi == 0) = [];
    Cx_idxi(Cx_idxi == 0) = [];
    Cx = sparse(Cx_idxi, Cx_idxj, Cx_vals, Nx*Nx, Nx*Nx);

    Cy_vals(Cy_idxi == 0) = [];
    Cy_idxj(Cy_idxi == 0) = [];
    Cy_idxi(Cy_idxi == 0) = [];
    Cy = sparse(Cy_idxi, Cy_idxj, Cy_vals, Nx*Nx, Nx*Nx);

    Dx_vals(Dx_idxi == 0) = [];
    Dx_idxj(Dx_idxi == 0) = [];
    Dx_idxi(Dx_idxi == 0) = [];
    Dx = sparse(Dx_idxi, Dx_idxj, Dx_vals, Nx*Nx, Nx*Nx);

    Dy_vals(Dy_idxi == 0) = [];
    Dy_idxj(Dy_idxi == 0) = [];
    Dy_idxi(Dy_idxi == 0) = [];
    Dy = sparse(Dy_idxi, Dy_idxj, Dy_vals, Nx*Nx, Nx*Nx);

    invEpsz = sparse(1:Nx*Nx, 1:Nx*Nx, 1./Epsz);
    Epsx = sparse(1:Nx*Nx, 1:Nx*Nx, Epsx);
    Epsy = sparse(1:Nx*Nx, 1:Nx*Nx, Epsy);

    Ax = Ax/dx; 
    Bx = Bx/dx;
    Cx = Cx/dx; 
    Dx = Dx/dx;
    Ay = Ay/dx; 
    By = By/dx;
    Cy = Cy/dx; 
    Dy = Dy/dx;

    %% Qxx, Qyy, Qxy, Qyx

    fprintf('Calculating Qs...\n');

    Qxx = -k0^(-2)*Ax*Dy*Cx*invEpsz*By + (Epsy + k0^(-2)*Ax*Dx)*(k0^2*I+Cy*invEpsz*By);
    Qyy = -k0^(-2)*Ay*Dx*Cy*invEpsz*Bx + (Epsx + k0^(-2)*Ay*Dy)*(k0^2*I+Cx*invEpsz*Bx);
    Qxy = k0^(-2)*Ax*Dy*(k0^2*I + Cx*invEpsz*Bx) - (Epsy + k0^(-2)*Ax*Dx)*Cy*invEpsz*Bx;
    Qyx = k0^(-2)*Ay*Dx*(k0^2*I + Cy*invEpsz*By) - (Epsx + k0^(-2)*Ay*Dy)*Cx*invEpsz*By;

    Q = [Qxx Qxy; Qyx Qyy];

    %% Diagonalisation

    fprintf('Taking Eigenvalues and Eigenvectors...\n');

    [eigvectors, eigvalues] = eigs(Q, NoModes, beta^2);
    beta = sqrt(diag(eigvalues));

    %% Ex, Ey, Ez

    fprintf('Calculating Ex, Ey, Ez, Hx, Hy, Hz...\n');

    Ex = zeros(Nx*Nx, NoModes);
    Ey = zeros(Nx*Nx, NoModes);
    Ez = zeros(Nx*Nx, NoModes);
    Hx = zeros(Nx*Nx, NoModes);
    Hy = zeros(Nx*Nx, NoModes);
    Hz = zeros(Nx*Nx, NoModes);
    for i = 1:NoModes
        Hx(:,i) = eigvectors(1:Nx*Nx, i);
        Hy(:,i) = eigvectors(Nx*Nx+1:2*Nx*Nx,i);
        Ez(:,i) = invEpsz*(-Dy*Hx(:,i) + Dx*Hy(:,i))/1i/w/eps0;
        Ey(:,i) = (-1i*w*mu0*Hx(:,i) - Ay*Ez(:,i))/1i/beta(i);
        Ex(:,i) = (1i*w*mu0*Hy(:,i) - Ax*Ez(:,i))/1i/beta(i);
        Hz(:,i) = -(-By*Ex(:,i) + Bx*Ey(:,i))/1i/w/mu0;    
    end

    %% Results

    for i = 1:NoModes
        RetVal.Ex{i} = ColumnToMatrix(Ex(:,i), Nx, Nx);
        RetVal.Ey{i} = ColumnToMatrix(Ey(:,i), Nx, Nx);
        RetVal.Ez{i} = ColumnToMatrix(Ez(:,i), Nx, Nx);
        RetVal.Hx{i} = ColumnToMatrix(Hx(:,i), Nx, Nx);
        RetVal.Hy{i} = ColumnToMatrix(Hy(:,i), Nx, Nx);
        RetVal.Hz{i} = ColumnToMatrix(Hz(:,i), Nx, Nx);
        RetVal.Eabs{i} = sqrt(abs(RetVal.Ex{i}).^2 + abs(RetVal.Ey{i}).^2 + abs(RetVal.Ez{i}).^2);
        RetVal.Habs{i} = sqrt(abs(RetVal.Hx{i}).^2 + abs(RetVal.Hy{i}).^2 + abs(RetVal.Hz{i}).^2);
    end
    RetVal.beta = beta;
    RetVal.nx = nx;
    RetVal.ny = ny;
    RetVal.nz = nz;
    RetVal.dx = dx;
    RetVal.lambda = lambda;
    RetVal.k0 = k0;
    RetVal.Nx = Nx;
    RetVal.PML_Depth = PML_Depth;
    RetVal.PML_TargetLoss = PML_TargetLoss;
    RetVal.PML_PolyDegree = PML_PolyDegree;
    RetVal.PML_SigmaMax = PML_SigmaMax;

end

%% Helper Functions

function M = ColumnToMatrix(C, Nx, Ny)

    M = reshape(C, Nx, Ny).';

end
function C = MatrixToColumn(M)
    
    M = M.';
    C = M(:);

end
