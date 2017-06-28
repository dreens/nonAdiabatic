%% Load E-field data for hop checking.
% nonAd.m loads data from nonadiabatic.dat generated by COMSOL
function hopgrid = nonAdClusterSwitching(sign)

%fprintf('Loading File: ')
all = importdata('nonadiabatic_switching.dat',' ',9);
all = all.data;

%fprintf('Reading xyz columns: ')
x = all(:,1);
y = all(:,2);
z = all(:,3);

% x,y have same range, z different. Get lookups that tell you, given
% an x vale x0, what ordinal it is, e.g. the 5th x-value.
ux = unique(x);
np = length(ux);
uz = unique(z);
npz = length(uz);
ntp = np*np*npz;
sp = max(diff(ux));
spz = max(diff(uz));
mp = min(ux);
mpz = min(uz);
l = @(x0) int16((x0 - mp)/sp + 1);
lz = @(z0) int16((z0 - mpz)/spz + 1);

% Get a single dimension matrix lookup index out of a triple index, for
% each row of the input data file.
fprintf('Loading linear indices: ')
xyz = sub2ind([np np npz],l(x),l(y),lz(z));

% Initialize all matrices
Ex = zeros([np np npz]);
Ey = zeros([np np npz]);
Ez = zeros([np np npz]);
xx = zeros([np np npz]);
yy = zeros([np np npz]);
zz = zeros([np np npz]);


% Fill all matrices
% fprintf('Filling Matrices: ')
xx(xyz) = x;
yy(xyz) = y;
zz(xyz) = z;
Ex(xyz) = all(:,4);
Ey(xyz) = all(:,5);
Ez(xyz) = all(:,6);

% Remove NaNs
nn(xyz) = isnan(Ex);
fprintf('NaNs: %d',sum(nn(:)))

Ex(nn) = 0;
Ey(nn) = 0;
Ez(nn) = 0;

% Save as first "p" point
Exp = Ex;
Eyp = Ey;
Ezp = Ez;
Enp = sqrt(Exp.^2+Eyp.^2+Ezp.^2);


%% Now do it again in the other switch state:

%fprintf('Loading File: ')
all = importdata('nonadiabatic_switching_s.dat',' ',9);
all = all.data;

%fprintf('Reading xyz columns: ')
x = all(:,1);
y = all(:,2);
z = all(:,3);

% x,y have same range, z different. Get lookups that tell you, given
% an x vale x0, what ordinal it is, e.g. the 5th x-value.
ux = unique(x);
np = length(ux);
uz = unique(z);
npz = length(uz);
ntp = np*np*npz;
sp = max(diff(ux));
spz = max(diff(uz));
mp = min(ux);
mpz = min(uz);
l = @(x0) int16((x0 - mp)/sp + 1);
lz = @(z0) int16((z0 - mpz)/spz + 1);

% Get a single dimension matrix lookup index out of a triple index, for
% each row of the input data file.
fprintf('Loading linear indices: ')
xyz = sub2ind([np np npz],l(x),l(y),lz(z));

% Initialize all matrices
Ex = zeros([np np npz]);
Ey = zeros([np np npz]);
Ez = zeros([np np npz]);
xx = zeros([np np npz]);
yy = zeros([np np npz]);
zz = zeros([np np npz]);


% Fill all matrices
% fprintf('Filling Matrices: ')
xx(xyz) = x;
yy(xyz) = y;
zz(xyz) = z;
Ex(xyz) = all(:,4);
Ey(xyz) = all(:,5);
Ez(xyz) = all(:,6);

% Remove NaNs
nn(xyz) = isnan(Ex);
fprintf('NaNs: %d',sum(nn(:)))

Ex(nn) = 0;
Ey(nn) = 0;
Ez(nn) = 0;

% Call these the "q" points
Exq = Ex*sign;
Eyq = Ey*sign;
Ezq = Ez*sign;
Enq = sqrt(Exq.^2+Eyq.^2+Ezq.^2);

%% Now we need to infer by symmetry...
% To get the fields after switching.
% 
% Ex2 = zeros([np np npz]);
% Ey2 = zeros([np np npz]);
% Ez2 = zeros([np np npz]);
% 
% Ex2(:,:,end:-1:1) = Ey;
% Ey2(:,:,end:-1:1) = -Ex;
% Ez2(:,:,end:-1:1) = -Ez;
% Ez2 = permute(Ez,[2 1 3]);
% 
% %typing convenience
% quivers = @(a,b,c,d) quiver(squeeze(a),squeeze(b),squeeze(c),squeeze(d));
% 
% % Okay that was some serious reflection. Let's confirm it worked.
% figure;
% subplot(1,2,1)
% p = 1;
% quivers(xx(:,p,:),zz(:,p,:),Ex(:,p,:),Ez(:,p,:))
% subplot(1,2,2)
% quivers(xx(:,p,:),zz(:,p,:),Ex2(:,p,:),Ez2(:,p,:))
%
% FORGET SYMMETRY. TOO ERROR PRONE. TOO HARD.
%

%% A test before full batch.
%hopgrid = zeros([np np npz]);
%k = sub2ind(size(hopgrid),np,np,10);
%checkHopSwitches(Exp(k),Eyp(k),Ezp(k),-Exq(k),-Eyq(k),-Ezq(k),100e-9);

%% Now we Loop through everything
hopgrid = zeros([np np npz]);
parfor k=1:np*np*npz
    hopgrid(k) = checkHopSwitches(Exp(k),Eyp(k),Ezp(k),Exq(k),Eyq(k),Ezq(k),100e-9);
end

end

%% Back from Cluster Switches
%
% Visualizing 3D grid of switch data as a fly through along z, indexed by
% phase angle.
%
[xxx yyy] = meshgrid(-1:0.1:1,-1:0.1:1);
absmax = max(Enp(:));
for i=1:2*npz-1
    figure;
    if i>npz
        j = 2*npz - i;
    else
        j = i;
    end
    mm = minus(:,:,j);
    pp = plus(:,:,j);
    full = [mm(end:-1:2,end:-1:2) pp(end:-1:2,:) ; ...
            pp(:,end:-1:2) mm];
    contourf(xxx,yyy,log10(full),'LineStyle','none')
    hold on
    exp = Exp(:,:,j);
    eyp = Eyp(:,:,j);
    exq = Exq(:,:,j);
    eyq = Eyq(:,:,j);
    fex = [exp(end:-1:2,end:-1:2), -exp(end:-1:2,:) ; ...
            -exp(:,end:-1:2), exp];
    fey = [eyp(end:-1:2,end:-1:2), eyp(end:-1:2,:) ; ...
            eyp(:,end:-1:2), eyp];

    fex2 = [exq(end:-1:2,end:-1:2), exq(end:-1:2,:) ; ...
            exq(:,end:-1:2), exq];
    fey2 = [eyq(end:-1:2,end:-1:2), -eyq(end:-1:2,:) ; ...
            -eyq(:,end:-1:2), eyq];
    hold on
    s = @(mat) mat(2:2:end,2:2:end);
    scale = 2/3 - log10(absmax/max(max(sqrt(fex.^2+fey.^2))))/4;
    q1 = quiver(s(xxx),s(yyy),s(-fex),s(-fey),scale);
    q1.Color = [0 .85 .2];
    q1.LineWidth = 3;
    scale = 2/3 - log10(absmax/max(max(sqrt(fex2.^2+fey2.^2))))/4;
    quiver(s(xxx),s(yyy),s(fex2),s(fey2),scale,'r','LineWidth',3)
    caxis([-5 -1]) 
    cmap = colormap(gray);
    colormap(flipud(cmap(20:end,:)))
    xlabel('X-axis (mm)','FontSize',13)
    ylabel('Y-axis (mm)','FontSize',13)
    phi = (i-1)*180/(npz-1)-90;
    title(sprintf('Switching Spin-Flips, \\phi=%3d^\\circ',phi),'FontSize',14)
    h = colorbar;
    ylabel(h,'log10 Hopping Probability','FontSize',13)
end
