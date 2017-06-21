%% Landau Zener Diabatic Transitions
% Using the 8x8 OH Hamiltonian.
% 
% Modified from TISE.m, an attempt to check Majorana loss in an octupole
% trap, see R:\COMSOL files\Magnetoelectrostatic CLoverleaf\
%
% Dave Reens, 11/24/15
%
%
%
% Modified again from LZ8state.m. In that file, the B-field is always
% oriented along the same axis. Therefore the rotation of the magnetic
% field and it's effect on the coupling between states is not treated. My
% suspicion is that this won't change the results too much, since the
% strongest coupling between states comes from the rotation of the relative
% angle between E and B, which is treated correctly here, but we'll see.
%
% Dave Reens, 1/18/17

h = 6.626e-34;
hb = h/(2*pi);
uB = 9.27401e-24;
uOH = 1.4*uB;
mOH = 17*1.67e-27;

% Say perfect quadrupole trap, homogeneous bias Efield, linear approach
% along line parallel to z-axis but offset by some minimum radius.
G = 400;        % Bfield Gradient,      T/m
v = 5;          % Molecule Velocity,    m/s
r = 25e-6;      % Minimum Radius,       m
E = 2.1e5;      % Electric Field,       V/m     (3kV/cm applied gives 2.1 in center)
IE = 2e7*h;     % Initial Energy,       J       (specifies how far from crossing to start)
ti = IE/uOH/G/v;% Initial time,         s

% Get Bfield, EB angle as functions of time to feed into Hamiltonian
Bz = @(t) -G*v*t;   
Br = @(t) G*r/2;
%Th = @(t) acos(Bz(t)./sqrt(Bz(t).^2+Br(t).^2)) ; % Theta goes smoothly from 0 to pi.
%H = @(t) OH_Ham_Simple_SI(sqrt(Bz(t).^2+Br(t)^2),E,Th(t));
H = @(t) OH_Ham_Lab_Fixed(Br(t),0,Bz(t),0,0,E);

% Set inputs to the ODE solver
t0 = 1e-8:1e-8:ti;
t0 = [fliplr(-t0) 0 t0];
[V,D] = eig(H(t0(1)));  % get initial eigenstates to initialize molecule
y0 = V(:,end);          % initialize molecule's state vector to highest energy eigenstate
y0 = y0 * 1i;           % Add phase to make sure it doesn't change solution

% This function modifies ode45's builtin plotting function to handle
% the complex valued state vector:
compmag = @(t,y,flag,varargin) odeplot(t,y.*conj(y),flag,varargin);

% The tolerance settings have been tuned for convergence. 1e-6, 1e-8 work.
odeop = odeset('RelTol',1e-6,'AbsTol',1e-8,'OutputFcn',compmag);

% ode45 takes a function that gives y'. For schrodinger eqn, y' = H*y/(i*hb)
tic
figure(1)
[ts, ys] = ode45(@(t,y) H(t)*y/hb/1i, t0, y0,odeop);
toc

% overwrite plot output during ode solving to remove annoying markers.
plot(ts,ys.*conj(ys),'LineWidth',2)

% Lets recast the state vector by independently diagonalizing the 
% hamiltonian at every time t
zs = zeros(size(ys));
es = zs;
for i=1:length(ts)
    [V,D] = eig(H(ts(i)));
    zs(i,:) = ys(i,:)*V;
    es(i,:) = diag(D);
end
close(figure(2))
figure(2);
cmap = get(gca,'ColorOrder');
set(gca,'ColorOrder',[cmap ; 0 0 0]); hold on
plot(ts,fliplr(zs.*conj(zs)),'LineWidth',2)
legend('Location','NorthEastOutside',...
    'f 3/2','f -3/2','f 1/2','f -1/2',...
    'e 1/2','e -1/2','e 3/2','e -3/2')
grid on
ylim([0 1])

%% Caught a significant bug.
% The Hamiltonian was changing basis at the LZ crossing, due to the angle
% changing sign from +pi/2 to -pi/2. We're used to ignoring this sign
% because it doesn't change the eigenvalues, but it does change the sign on
% some of the eigenvectors. Such a sudden basis change is unacceptable,
% unless we were to also change the basis of the statevector y together
% with the hamiltonian during the ODE solver. I fixed this by adding an
% absolute value on Th(t).
%
% Actually, adding absolute value still creates an unphysical kink. I
% changed from atan to acos, so that theta goes nicely from 0 to pi,
% smoothly passing pi/2 at the loss plane.
%
% Notably, there are now no longer hops to any state other than f -3/2.

%% Compare to formula
fprintf('Schrodinger Probability is %2.1f%%\n',100*abs(zs(end,7))^2)

middle = (length(ts)+1)/2;
gap = es(middle,8)-es(middle,7);
fprintf('Gap width is %1.2f MHz\n',gap/h/1e6)
dhdt = 400*uOH*v*2;
ex = pi*gap^2/(2*dhdt*hb);
p = exp(-ex);
fprintf('Calculated Probability is %2.1f%%\n',100*p)



