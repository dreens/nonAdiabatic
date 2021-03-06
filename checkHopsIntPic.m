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

function out = checkHopsIntPic(z,Exp,Eyp,Ezp,v)

z   = squeeze(z);
Exp = squeeze(Exp);
Eyp = squeeze(Eyp);
Ezp = squeeze(Ezp);

h = 6.626e-34;
hb = h/(2*pi);
uB = 9.27401e-24;
uOH = 1.4*uB;
mOH = 17*1.67e-27;

Ex = griddedInterpolant(z,Exp,'cubic');
Ey = griddedInterpolant(z,Eyp,'cubic');
Ez = griddedInterpolant(z,Ezp,'cubic');

H0 = OH_Ham_Lab_Fixed(1e-4,0,0,Ex(min(z)),Ey(min(z)),Ez(min(z)));
VS = @(t) (OH_Ham_Lab_Fixed(1e-4,0,0,Ex(v*t+min(z)),Ey(v*t+min(z)),Ez(v*t+min(z))) - H0);
VI = @(t) expm(1i*H0*t/hb)*VS(t)*expm(-1i*H0*t/hb);


% Set inputs to the ODE solver
t0 = 0:1e-9:((max(z)-min(z))/v);
[V,D] = eig(H0);  % get initial eigenstates to initialize molecule
y0 = V(:,end);          % initialize molecule's state vector to highest energy eigenstate
y0 = y0 * 1i;           % Add phase to make sure it doesn't change solution

% This function modifies ode45's builtin plotting function to handle
% the complex valued state vector:
compmag = @(t,y,flag,varargin) odeplot(t,y.*conj(y),flag,varargin);

% The tolerance settings have been tuned for convergence. 1e-6, 1e-8 work.
odeop = odeset('RelTol',1e-8,'AbsTol',1e-10,'OutputFcn',compmag);

% ode45 takes a function that gives y'. For schrodinger eqn, y' = H*y/(i*hb)
figure(1)
[ts, ys] = ode45(@(t,y) VI(t)*y/1i/hb, t0, y0, odeop);

% overwrite plot output during ode solving to remove annoying markers.
plot(ts,ys.*conj(ys),'LineWidth',2)

% Lets recast the state vector by independently diagonalizing the 
% hamiltonian at every time t
zs = zeros(size(ys));
xxs = zs;
es = zs;
for i=1:length(ts)
    xxs(i,:) = (expm(-1i*H0*ts(i)/hb)*(ys(i,:)'))';
    [V,D] = eig(H0+VS(ts(i)));
    zs(i,:) = conj(xxs(i,:))*V;
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

out = ys;

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


