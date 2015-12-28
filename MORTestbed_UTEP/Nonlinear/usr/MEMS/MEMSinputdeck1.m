function  [xold,consts,operators] = MEMSinputdeck1(nodes,DLim)
%This function sets up the MEMS function (constants, operations, error
%function, plot function).
%--------------------------------------------------------------------------
%Inputs:
%-------
%N         - scalar indicating the problem discretization in the
%            x-direction
%M         - scalar indicating the problem discretization in the
%            y-direction
%
%Outputs:
%--------
%xold      - N*(M+2) x 1 solution ("guess" or current iterate) vector 
%consts    - vector containing the relevant constants of the problem
%            (the format is [S0 E In0 rho0 lambda p0 z0 width dy,...
%            viscosity, dx,height b0 dt M N]
%operators - cell array containing the relevant operators for the
%            problem (operators = [{INT} {E} {D2} {grady} {B} {C},...
%            {gradx2} {gradx} {LAP} {D4}]) 
%--------------------------------------------------------------------------

%Load all constant parameters for the system
%consts = [S0 E In0 rho0 lambda p0 z0 width dy viscosity dx height b0 dt];
M = nodes(1);
N = nodes(2);

%Physical geometry parameters (see Reweinski Thesis)
del = 1e-7;
height = 2.2; 
long = DLim(1,2) - DLim(1,1);
width = DLim(2,2) - DLim(2,1);
z0=2.3;
dx = (long)/(N+1);
dy = (width)/(M+1);
% dx=long/(N+1); 
% dy=width/(M+1);

amp = 1;

%Material and air parameters (see Reweinski Thesis)
epsilon0=8.854e-6;
% rho0=2330e-18;
rho0=2300e-18;
E=149.0e3*amp;
S0=-3.7*amp;
In0=1/12; 
atm_Pa_conversion=1.013e-1;
viscosity=1.82e-11;
lambda=0.064;
p0=atm_Pa_conversion;
% p0 =  1.0885764*atm_Pa_conversion;

% Input parameter (see Reweinski Thesis)
b0 = -(3*epsilon0*del/(2*rho0*height));

%Store these constants in a single vector for transfer between functions
consts = [S0 E In0 rho0 lambda p0 z0 width dy viscosity dx height b0 del];

%Store the initial condition for the MEMS problem
xold = zeros(N*(M+2),1); xold(1:N) = z0; xold((2*N+1):((M+2)*N)) = p0;

%Set up operators
% Operator INT
IntSeg = ones(1,M);
%Boundary points are p=0 here
IntSeg = dy*IntSeg;
IntRow = zeros(1,M*N); IntRow(1:M) = IntSeg;
for i = 1:N
    INT(i,:) = IntRow;
    IntRow = circshift(IntRow,[0 M]);
end
INT=sparse(INT); 

%Operator D4
D4 = 6*diag(ones(N,1))-4*diag(ones(N-1,1),1)-4*diag(ones(N-1,1),-1)+diag(ones(N-2,1), 2)+diag(ones(N-2,1),-2);
D4(1,1)=D4(1,1)+1;
D4(N,N)=D4(N,N)+1;
D4=(1.0/dx^4)*D4;
D4 = sparse(D4);

%Operator D2
D2 = -2*diag(ones(N,1))+diag(ones(N-1,1),1)+diag(ones(N-1,1),-1);
D2(1,1)=-1;
D2(N,N)=-1;
D2=(1.0/dx^2)*D2;
D2 = sparse(D2);

%grady
grady = zeros(M*N,M*N);
grady(1:M,1:M) = (1.0/(2.0*dy))*(diag(ones(M-1,1),1)-diag(ones(M-1,1), -1));
for index=2:N
    grady([((index-1)*M+1):(index*M)], [((index-1)*M+1):(index*M)]) = ...
        grady([1:M],[1:M]);
end;
grady=sparse(grady);

% gradx
gradx = zeros(M*N,M*N);
e = ones((N+1)*M,1)/(2.0*dx);
%   gradx = diag(e, M) - diag(e, -M);
gradx = spdiags([e -1.0*e],[M -M],(M*N),(M*N));
% -spdiags(ones((N-1)*M,1), -M,(M*N)*(M*N)));
for index = 1:M
    gradx(index,index) = -1.0/(2*dx);
    gradx(M*(N-1)+index, M*(N-1)+index) = 1.0/(2*dx);
end;
gradx = sparse(gradx);

%Laplacian
e=1/dx^2; f1=-2*1/dx^2; f2=-2*1/dy^2; f=f1+f2; g=1/dy^2;
diag0=ones(N*M,1);
diag1=ones(N*M-1,1);
%Boundary cond at 0 & L, dp/dx=0 => p_1=p_0 => d^2/dx^2 p_1= 1/dx^2 * (-x1 + x2)
diag0(1:M)=1/2;
diag0((N-1)*M+1:N*M)=1/2;
for i=1:N-1
    diag1(i*M)=0;
end
LAP=f1*diag(diag0)+f2*diag(ones(N*M,1))+g*diag(diag1,1)+g*diag(diag1,-1)+e*diag(ones((N-1)*M,1),M)+e*diag(ones((N-1)*M,1),-M);
LAP=sparse(LAP);

%gradx2
gradx2 = (1/(2*dx))*(diag(ones(N-1,1),1)-diag(ones(N-1,1),-1));
gradx2=sparse(gradx2);

% Operator B -- input vector
B = zeros((M+2)*N,1);
B((N+1):(2*N)) = 1;
B = B*b0;
B=sparse(B);

% Operator C -- Output vector
C1=zeros(1,(M+2)*N); C2 = C1; C3 = C1;
C1(round(N/2))=1;   % Vertical displacement at center
C3(2*N+round(M*N/2)) = 1; % Pressure under beam center
C2(N+round(N/2)) = 1; % du^3/dt at center of beam
%C = [C1;C2;C3];
C = [C1];
C=sparse(C);
E = speye(length(B));

%Store all operators into cell array for transfer between functions
operators = [{INT} {E} {D2} {grady} {B} {C} {gradx2} {gradx} {LAP} {D4} {xold}];

end