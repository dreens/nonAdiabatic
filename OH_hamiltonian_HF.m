    %% used for OH ground hyperfine state calculation in E, B fields

function HFS = OH_hamiltonian_HF(Bf,E,beta)

%E = E; %input in V/m
sb = sin(beta);
cb = cos(beta);

uB=9.27401*1e-24; %this way input can be in Tesla
g = 0.9355; %1.4032 is 3/2 * the gJ factor for OH, 0.9355
uE=3.33564*1e-30; %so uE (C m) * E (V/m) gives Joules.
dE=1.67;
h=6.62607*1e-34;
c=2.99792458*1e8;
LD = h*(1.667358e9);

%% now the hyperfine parameters (in Hz)
Gamma = h*54.1e6;
Xi = h*1.96e6;
Gamma_e = Gamma - Xi/2;
Gamma_f = Gamma + Xi/2;


deltag=1.267*1e-3;
gf=g-deltag/2;
ge=g+deltag/2;

beta_e=uB*ge*Bf;
beta_f=uB*gf*Bf;

%This 2x2 gives the coupling between F=2 and F=1 hyperfine states.
A = @(ga,b,mF) [ga/2 + 2*b*mF*3/8      b*sqrt(4-mF^2)/4        ;
               b*sqrt(4-mF^2)/4      -ga/2 + 2*b*mF*5/8       ];

%Now we assemble these for each mF into 2 8x8's for e,f parity
Af = @(ga,b) [ ga/2+b*3/2    zeros(1,7)                                           ;
               zeros(2,1)    A(ga,b,1)   zeros(2,5)                               ;
                             zeros(2,3)  A(ga,b,0)    zeros(2,3)                  ;
                                         zeros(2,5)   A(ga,b,-1)    zeros(2,1)    ;
                                                      zeros(1,7)    ga/2-b*3/2   ];
 

%And finally into a 16x16 including parity as well
F = [Af(Gamma_f,beta_f)            zeros(8)      ;
          zeros(8)            Af(Gamma_e,beta_e)  ];

%% Now we add lambda doubling:
% Since LD is the known Lambda Doubling between the stretched states, we
% need to take away the Xi factor already splitting the energies via the
% hyperfine Gamma parameter.

F = F - (LD-Xi/2)*[zeros(8) zeros(8) ; zeros(8) eye(8)] - Gamma_f/2 * eye(16);
      
% here the energy zero is taken to be the F2 states @ B=0, hence the -Gamma_f/2.
      
%% Next we need to change basis into I-J

% Right now the matrix is in this order:
% 2 2
% 2 1
% 1 1
% 2 0
% 1 0
% 2 -1
% 1 -1
% 2 -2

% Encoding a CG matrix that takes IJ to F:
CG = zeros(8);
CG(1,1) = 1;
CG(2:3,2:3) = [1/4 3/4 ; 3/4 -1/4];
CG(4:5,4:5) = [1/2 1/2 ; 1/2 -1/2];
CG(6:7,6:7) = [3/4 1/4 ; 1/4 -3/4];
CG(8,8) = 1;
CG = arrayfun(@(x) sqrt(x)*(x>0) - sqrt(-x)*(x<0),CG);

CG16 = [CG zeros(8) ; zeros(8) CG];

FIJ = CG16 * F * CG16;


%% Finally we add in the stark hamiltonian
%A = eye(4).*g*repmat([-1.5 -.5 .5 1.5]',1,4);
%A = A*uB*Bf - LD*eye(4)/2;
%D = A + LD*eye(4);
B = eye(4).*repmat([0.6 .2 -.2 -0.6]',1,4)*uE*E*cb*dE;
BB = [0 3 0 0 ; 3 0 4 0 ; 0 4 0 3 ; 0 0 3 0];
BB = -0.2*sqrt(BB)*uE*E*sb*dE;
B = B + BB;
B = rot90(B,2);
S  = [zeros(4) B ; B zeros(4)];
SH = [S zeros(8) ; zeros(8) S];
SH = SH(:,[1 9 2 10 3 11 4 12 5 13 6 14 7 15 8 16]);
SH = SH([1 9 2 10 3 11 4 12 5 13 6 14 7 15 8 16],:);

HFS = FIJ+SH;




