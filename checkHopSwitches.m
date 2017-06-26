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
%
% Modified again for purposes of studying nonadiabatic transitions in the
% decelerator.

function hop = checkHopSwitches(Exp,Eyp,Ezp,Exq,Eyq,Ezq,tau)

h = 6.626e-34;
hb = h/(2*pi);
uB = 9.27401e-24;
uOH = 1.4*uB;
mOH = 17*1.67e-27;
a = @(t) exp(-t/tau);
b = @(t) 1-a(t);
H = @(t) OH_Ham_Lab_Fixed(1e-10,0,0,Exp*a(t)+Exq*b(t),Eyp*a(t)+Eyq*b(t),Ezp*a(t)+Ezq*b(t));

% Set inputs to the ODE solver
t0 = 0:1e-9:3*tau;
[V,D] = eig(H(t0(1)));  % get initial eigenstates to initialize molecule
[~,l] = sort(diag(real(D)));
y0 = V(:,l(end));          % initialize molecule's state vector to highest energy eigenstate
y0 = y0 * 1i;           % Add phase to make sure it doesn't change solution

% This function modifies ode45's builtin plotting function to handle
% the complex valued state vector:
compmag = @(t,y,flag,varargin) odeplot(t,y.*conj(y),flag,varargin);

% The tolerance settings have been tuned for convergence. 1e-6, 1e-8 work.
odeop = odeset('RelTol',1e-8,'AbsTol',1e-10,'OutputFcn',compmag);

% ode45 takes a function that gives y'. For schrodinger eqn, y' = H*y/(i*hb)
figure(1)
[ts, ys] = ode45(@(t,y) H(t)*y/hb/1i, t0, y0, odeop);

% overwrite plot output during ode solving to remove annoying markers.
plot(ts,ys.*conj(ys),'LineWidth',2)

% Lets recast the state vector by independently diagonalizing the 
% hamiltonian at every time t
zs = zeros(size(ys));
es = zs;
for i=1:length(ts)
    [V,D] = eig(H(ts(i)));
    [Ds, l] = sort(real(diag(D)));
    temp = conj(ys(i,:))*V;
    zs(i,:) = temp(l);
    es(i,:) = real(Ds);
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

%full = ys;
final = zs.*conj(zs);
hop = sum(final(end,1:end-2));

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
% fprintf('Schrodinger Probability is %2.1f%%\n',100*abs(zs(end,7))^2)
% 
% middle = (length(ts)+1)/2;
% gap = es(middle,8)-es(middle,7);
% fprintf('Gap width is %1.2f MHz\n',gap/h/1e6)
% dhdt = 400*uOH*v*2;
% ex = pi*gap^2/(2*dhdt*hb);
% p = exp(-ex);
% fprintf('Calculated Probability is %2.1f%%\n',100*p)
% 


