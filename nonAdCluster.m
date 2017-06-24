%% Load E-field data for hop checking.
% nonAd.m loads data from nonadiabatic.dat generated by COMSOL
function hopgrid = nonAdCluster()

%fprintf('Loading File: ')
all = importdata('nonadiabatic.dat',' ',9);
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


%% Now we Loop through and check hopping probability
hopgrid = zeros(np);
zzz = zz(1,1,:);
parfor k=1:np*np
    i = mod(k-1,np) + 1;
    j = floor((k-1)/np)+1;
    hopgrid(k) = checkHops(zzz,Ex(i,j,:),Ey(i,j,:),Ez(i,j,:),750000);
end

end
