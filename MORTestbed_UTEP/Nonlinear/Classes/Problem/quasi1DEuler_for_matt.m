function [U] = quasi1DEuler()

% !! Only for an equispaced mesh 
% Physical constants
gamma = 1.4;
R = 287.06;
Cv = R / (gamma-1);
% Boundary Conditions
M0 = 0.5;
T0 = 300;
P0 = 10^6;
c0 = sqrt(gamma*R*T0);
rho0 = P0 / (R*T0);
u0 = M0 * c0;
Tt = T0*(1 + (gamma-1)/2*M0^2); % total temperature
Pt = P0 *(1 + (gamma-1)/2*M0^2)^(gamma/(gamma-1));% total pressure
% Tt = 300;
% Pt = 10^8; %in N/m^2 = Pa or 10^5
%nPoints = 41;
nPoints = 501;
%nPoints = 1001;
nVol = nPoints +1;

% Create the mesh
xmax = 2; % to be changed
p = 0.4;
[x,y] = inletMesh(nPoints,xmax,p); % function to be changed
dx = x(2,1) - x(1,1);
% Compute S(x) and SVol(x)
[S,SVol,xVol] = computeNozzleArea(x,y,nVol);

% Compute dt
 dt = 8*dx / (u0 + c0);

% Initial Conditions
M = M0*ones(nVol,1);
rho = rho0*ones(nVol,1);
P = P0*ones(nVol,1);
u = u0*ones(nVol,1);
% plots
figure(1);
hold on
plot(xVol,M,'-o')
axis([-0.0125 xmax*1.0125 0 2.5])
figure(2)
hold on
plot(xVol,P/P0,'r')
axis([-0.0125 xmax*1.0125 0 2.5])

% Build the state vector U
U = primitiveToConservative(rho,u,P,gamma);
iIter = 0;
dU = 0;
tol = 10^(-6);
% Iterations
while(norm(dU)>tol*norm(U) || iIter==0) % termination criterion to be changed
    iIter = iIter + 1;
    fprintf('\n       **** Iteration = %g  **** \n\n',iIter);
    
    % Boundary Conditions
    [Ainlet,Binlet,dVinlet] = semiImplicitInletBC(U,SVol,dt,dx,gamma,Cv,Pt,Tt);
    % at the outlet
    dVoutlet = semiImplicitOutletBC(U,SVol,dt,dx,gamma,Cv);
    % Build the Fluxes
    F = fluxes(U,gamma);
    % Build Roe matrices
    [Jac] = roeMat(U,gamma);
    % Build Roe fluxes
    roeF = roeFlux(U,Jac,F);
    %Extract the pressure from U
    P = extractP(U,gamma);
    % Build the RHS
    RHS = buildRHS(dVinlet,dVoutlet,dx,dt,roeF,S,SVol,P);
    % Build the LHS
    LHS = buildLHS(dx,dt,gamma,S,SVol,U,Jac,Ainlet,Binlet);
    % solve for the increment
    dU = matrixSolve(LHS,RHS);
    % Update the state vector
    U = U + dU;
    if (mod(iIter,40) == 0)
        [rho,u,P,M] = conservativeToPrimitive(U,gamma);
        figure(2)
        plot(xVol,P/P0,'r-o');
        figure(1)
        plot(xVol,M,'r-o');
    end
end
figure(1)
plot(xVol,M,'k-o');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [S,SVol,xVol] = computeNozzleArea(x,y,nVol)
    nPoints = nVol-1;
    SVol = zeros(nVol,1);
    xVol = zeros(nVol,1);
    S = abs(y);
    SVol(1,1) = S(1,1);
    SVol(nVol,1) = S(nPoints,1);
    SVol(2:(nVol)-1,1) = 0.5*(S(1:nVol-2,1) +S(2:(nVol-1),1));
    xVol(1:(nVol-1),1) = x(1:(nVol-1),1) - (x(2,1)-x(1,1))/2*ones(nVol-1,1);
    xVol(nVol,1) = x(nVol-1,1) + (x(2,1) -x(1,1))/2;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function U = primitiveToConservative(rho,u,P,gamma)

U(:,1) = rho;
U(:,2) = rho.*u;
U(:,3) = P/(gamma-1) + rho.*u.^2/2;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [rho,u,P,M] = conservativeToPrimitive(U,gamma)
      rho = U(:,1);
      u = U(:,2)./U(:,1);
      P = (gamma-1) * (U(:,3) - rho.*u.^2/2);
      M = u./ (sqrt(gamma*P./rho));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,A] = inletMesh(nPoints,xmax,p)

   x = linspace(0,xmax,nPoints)';
   A = (0.6*(x-1).^2 + p);


figure(3)
hold on
plot(x,A)
for i=1:nPoints
    plot([x(i,1) x(i,1)],[0 A(i,1)]);
end
plot([0,xmax],[0,0]);
axis([0 xmax 0 max(A)])

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dV = semiImplicitOutletBC(U,SVol,dt,dx,gamma,Cv)

[nPoints,dim] = size(U);
dV = zeros(dim,1);
[rhoO,uO,PO,MO] = conservativeToPrimitive(U(nPoints-1:nPoints,:),gamma);
cO = uO./MO;
SO(1,1) = SVol(nPoints-1,1);
SO(2,1) = SVol(nPoints,1);
% average values
rho = 0.5*(rhoO(1,1) + rhoO(2,1));
u = 0.5*(uO(1,1) + uO(2,1));
P = 0.5*(PO(1,1) + PO(2,1));
c = 0.5*(cO(1,1) + cO(2,1));
% eigenvalues
lambda1 = u/dx;
lambda2 = (u+c)/dx;
lambda4 = (u-c)/dx;
R1 = -lambda1 /(lambda1+1/dt) * (rhoO(2,1) - rhoO(1,1) -1/c^2*(PO(2,1)-PO(1,1)));
R2 = -lambda2 / (lambda2+1/dt) * (PO(2,1) - PO(1,1) + rho*c*(uO(2,1)-uO(1,1))) - gamma*P*u / (1/dt+lambda2) * (SO(2,1)-SO(1,1)) /(SO(2,1)*dx);
R4 = -lambda4 / (lambda4+1/dt) * (PO(2,1) - PO(1,1) - rho*c*(uO(2,1)-uO(1,1))) - gamma*P*u / (1/dt+lambda4) * (SO(2,1)-SO(1,1)) /(SO(2,1)*dx);
if (MO(1,1) > 1)
    dP = 0.5*(R2 + R4); % supersonic outlet
else
    dP = 0;
end
drho = R1 + dP / c^2;
du = (R2 - dP) / (rho*c);
T = (PO(2,1) + dP) / ( (gamma-1)*Cv*(rhoO(2,1) + drho));
e = (rhoO(2,1)+drho) * ( Cv*T + 0.5*(uO(2,1) + du)^2);
dV(1,1) = rhoO(2,1) + drho - U(nPoints,1);
dV(2,1) = (rhoO(2,1)+drho)*(uO(2,1)+du) - U(nPoints,2);
dV(3,1) = e - U(nPoints,3);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,B,dV] = semiImplicitInletBC(U,SVol,dt,dx,gamma,Cv,Pt,Tt)

[nPoints,dim] = size(U);
[rhoI,uI,PI,MI] = conservativeToPrimitive(U(1:2,:),gamma);
if (MI(1,1) > 1)
    A = eye(dim);
    B = zeros(dim,dim);
    dV = zeros(dim,1);
else
    cI = uI./MI;
    % average values
    MI(1,1)
    rho = 0.5*(rhoI(1,1) + rhoI(2,1));
    u = 0.5*(uI(1,1) + uI(2,1));
    c = 0.5*(cI(1,1) + cI(2,1));
    cstar = sqrt(2*gamma*(gamma-1)*Cv*Tt/(gamma+1));
    dPdu = -2*gamma*uI(1,1)*Pt/((gamma+1)*cstar^2)*(1-(gamma-1)/(gamma+1)*(uI(1,1)^2)/(cstar^2))^(1/(gamma-1));
    dTdu = -2*Tt*(gamma-1)/((gamma+1)*cstar^2)*uI(1,1);
    TI = Tt*(1-(gamma-1)/(gamma+1)*(uI(1,1)^2)/cstar^2);
    
    drhodu = 1/((gamma-1)*Cv*TI^2)*(TI*dPdu - PI(1,1)*dTdu);
    lambda3 = (u-c)/dx;
    R3 = -(u-c)/dx*(PI(2,1)-PI(1,1)-rho*c*(uI(2,1)-uI(1,1))) - gamma*PI(1,1)*uI(1,1)*(SVol(2,1)-SVol(1,1))/(SVol(1,1)*dx);
    dV = [0;0;R3];
    a12 = -drhodu;
    a22 = - dPdu;
    a32 = -rho*c*(1/dt-(u-c)/dx) + gamma*PI(1,1)*(SVol(2,1)-SVol(1,1))/(SVol(1,1)*dx);
    a33 = 1/dt-(u-c)/dx + gamma*uI(1,1)*(SVol(2,1)-SVol(1,1))/(SVol(1,1)*dx);
    Aprime = [1 a12 0; ...
            0 a22 1; ...
            0 a32 a33];
    beta = gamma-1;
    alpha = u^2/2;

    iS = [1 0 0; u, rho 0; alpha rho*u 1/beta];
    if (u==0)
        sol = [0; ...
                  R3/(-rho*c*(1/dt-lambda3)+gamma*PI(1,1)*(SVol(2,1)-SVol(1,1))/(SVol(1,1)*dx));...
                  0];  
    else      
        sol = Aprime\dV;
    end
    dV = iS*sol;
    A = eye(dim);
    B = zeros(dim,dim);
end

end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dV] = matrixSolve(Mat,dU)

  [nVol,dim] = size(dU);
  rhs = zeros(nVol*dim,1);
  rhs(:) = dU';
  res = Mat\rhs; 
  
  dV = reshape(res,dim,nVol)';
 


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F = fluxes(U,gamma)

F = zeros(size(U));
F(:,1) = U(:,2);
F(:,2) = U(:,2).^2 ./ U(:,1) + (gamma-1)*(U(:,3)- 0.5 * U(:,2).^2 ./ U(:,1));
F(:,3) = U(:,2) ./ U(:,1) .* (gamma*U(:,3) - 0.5*(gamma-1)* U(:,2).^2 ./ U(:,1));

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Jac] = roeMat(U,gamma)

[nVol,~] = size(U);
Jac(nVol-1).Ap = 0;
beta = gamma-1;
P = (gamma-1) * (U(:,3) - 0.5* U(:,2).^2 ./ U(:,1));
V = 1:nVol-1;
rhoR = sqrt(U(V,1).*U(V+1,1));
uR = (U(V,2)./sqrt(U(V,1)) + U(V+1,2)./sqrt(U(V+1,1))) ./ (sqrt(U(V,1)) + sqrt(U(V+1,1)));
hR = ((U(V,3)+P(V,1))./sqrt(U(V,1)) + (U(V+1,3)+P(V+1,1))./sqrt(U(V+1,1))) ./ (sqrt(U(V,1)) + sqrt(U(V+1,1)));
cR = sqrt((gamma-1) * (hR -0.5*uR.^2));
[rho,u,P,M] = conservativeToPrimitive(U,gamma);
c = u./M;
lambda = [u,u+c,u-c];
sigma0 = 10; % constant for Roe entropy correction
% For every point
for i=1:nVol-1
    % Build the matrices S and CA
    alpha = uR(i,1)^2 / 2;
    S = [1 0 0; -uR(i,1)/rhoR(i,1) 1/rhoR(i,1) 0; alpha*beta -uR(i,1)*beta beta];
    iS = [1 0 0; uR(i,1) rhoR(i,1) 0; alpha rhoR(i,1)*uR(i,1) 1/beta];
    CA = [1 0 -1/cR(i,1)^2; 0 rhoR(i,1)*cR(i,1) 1; 0 -rhoR(i,1)*cR(i,1) 1];
    iCA = [1 1/(2*cR(i,1)^2) 1/(2*cR(i,1)^2); 0 1/(2*rhoR(i,1)*cR(i,1)) -1/(2*rhoR(i,1)*cR(i,1)); 0 1/2 1/2];    
    lambdaR = [uR(i,1) 0 0; 0 uR(i,1)+cR(i,1) 0; 0 0 uR(i,1)-cR(i,1)];
    for j=1:3  %  entropy correction
        eps = sigma0*max([0,lambdaR(j,j)-lambda(i,j),lambda(i+1,j)-lambdaR(j,j)]);
        if (abs(lambdaR(j,j))<eps)
            lambdaR(j,j) = 0.5*(lambdaR(j,j)^2/eps+eps);
        end
    end
    lambdaRp = (lambdaR + abs(lambdaR)) / 2;    
    A = iS*iCA*lambdaR*CA*S;
    Jac(i).Ap = iS*iCA*lambdaRp*CA*S;
    Jac(i).Am = A - Jac(i).Ap;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function roeF = roeFlux(U,Jac,F)

[nVol,dim] = size(U);
roeF =zeros(nVol-1,dim);
for i=1:nVol-1
    absA = Jac(i).Ap - Jac(i).Am;
    roeF(i,:) = (F(i,:) + F(i+1,:))/2 - (absA * (U(i+1,:) - U(i,:))' / 2)';
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = extractP(U,gamma)

P = (gamma-1) * (U(:,3) - 0.5* U(:,2).^2 ./ U(:,1));

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function RHS = buildRHS(dVinlet,dVoutlet,dx,dt,roeF,S,SVol,P)

dim = size(dVinlet,1);
nVol = size(roeF,1)+1;
RHS = zeros(nVol,dim);

% Inlet BC
RHS(1,1:dim) = dVinlet;

% Outlet BC
RHS(nVol,1:dim) = dVoutlet;
% Other blocks
V = 2:(nVol-1);
Svec = repmat(S(V,1),1,dim);
Svecm = repmat(S(V-1,1),1,dim);
SVolvec = repmat(SVol(V,1),1,dim);
RHS(V,1:dim) = -(dt ./  SVolvec) .* (Svec.*roeF(V,1:dim) - Svecm.*roeF(V-1,1:dim)) / dx + dt * [zeros(nVol-2,1), (S(V,1) - S(V-1,1)) .* P(V,1) ./ (SVol(V,1) * dx),zeros(nVol-2,1)];
end
%%%%%%%%%%%%%%%%%%%%
function [Mat] = buildLHS(dx,dt,gamma,S,SVol,U,Jac,Ainlet,Binlet)

[nVol,dim] = size(U);
Mat = sparse(dim*nVol,dim*nVol);
%block(nVol).A = 0;

% Inlet Boundary
centBlock = 1:dim;
nextBlock = (dim+1):2*dim;
Mat(centBlock,centBlock) = Ainlet;
Mat(centBlock,nextBlock) = Binlet;
%block(1).A = eye(dim);
%block(1).B = zeros(dim,dim);
%block(1).C = zeros(dim,dim);
% Outlet Boundary
centBlock = ((nVol-1)*dim+1):(nVol*dim);
Mat(centBlock,centBlock) = eye(dim);
%block(nVol).A = eye(dim);
%block(nVol).B = zeros(dim,dim);
%block(nVol).C = zeros(dim,dim);
for i=2:nVol-1
    Api = Jac(i).Ap;
    Ami = Jac(i).Am;
    Apim = Jac(i-1).Ap;
    Amim = Jac(i-1).Am;   
    
    H = (S(i,1) - S(i-1,1)) / (SVol(i,1)*dx) * (gamma-1) * [0 0 0; U(i,2)^2 / (2*U(i,1)^2) -U(i,2)/U(i,1) 1; 0 0 0];
    centBlock = ((i-1)*dim+1):(i*dim); 
    %Mat(centBlock,centBlock) = eye(dim) +  dt / (SVol(i,1) * dx) * (S(i,1) * Api - S(i-1,1)*Amim) -dt * H;
    Mat(centBlock,centBlock) = eye(dim) +  dt / (SVol(i,1) * dx) * (S(i,1) * Apim - S(i-1,1)*Ami) -dt * H;
    if (i<nVol)
        nextBlock = (i*dim+1):((i+1)*dim);
        Mat(centBlock,nextBlock) = dt / (SVol(i,1) * dx) * S(i,1) * Ami;
    end
    if (i>1)
        prevBlock = ((i-2)*dim+1):((i-1)*dim);
        Mat(centBlock,prevBlock) = - dt / (SVol(i,1) * dx) * S(i-1,1) * Apim; 
    end   
end

end