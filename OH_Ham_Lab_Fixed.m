%% OH ground state in E, B fields without fixing B in z direction
%
% I realized that our insistence on fixing the basis to have magnetic field
% in the z-direction is invalid when solving the time-independent
% Schrodinger equation for Landau-Zener hopping.
%
% This modified function is based on Lara, Lev, and Bohn's PRA 78.033433
% appendix A. It's just molecule case (a) rotational theory but they've
% pre-evaluated the Wigner-rotation terms into 3j symbols so it's a great
% startpoint.
%

% This function is working as of 9:20AM, 1/18/17. I have tested that equal
% magnitude E fields in different directions give the same eigenvalues.
% Same for B fields. I have confirmed that results are equivalent to those
% for OH_Ham_Simple_SI(B,E,theta).

function H = OH_Ham_Lab_Fixed(Bx,By,Bz,Ex,Ey,Ez)

persistent cacheH;

if isempty(cacheH)
    cacheH.full = 0;
    cacheH.H0 = OH_Ham_Lab_Fixed(0,0,0,0,0,0);
    cacheH.Bx = OH_Ham_Lab_Fixed(1,0,0,0,0,0) - cacheH.H0;
    cacheH.By = OH_Ham_Lab_Fixed(0,1,0,0,0,0) - cacheH.H0;
    cacheH.Bz = OH_Ham_Lab_Fixed(0,0,1,0,0,0) - cacheH.H0;
    cacheH.Ex = OH_Ham_Lab_Fixed(0,0,0,1,0,0) - cacheH.H0;
    cacheH.Ey = OH_Ham_Lab_Fixed(0,0,0,0,1,0) - cacheH.H0;
    cacheH.Ez = OH_Ham_Lab_Fixed(0,0,0,0,0,1) - cacheH.H0;
    cacheH.full = 1;
    H = OH_Ham_Lab_Fixed(Bx,By,Bz,Ex,Ey,Ez);
elseif cacheH.full
    H = cacheH.H0 + ...
        cacheH.Bx*Bx + cacheH.By*By + cacheH.Bz*Bz + ...
         cacheH.Ex*Ex + cacheH.Ey*Ey + cacheH.Ez*Ez;
else

    uB=9.27401*1e-24; 
    g = 0.9355; %1.4032 is 3/2 * the gJ factor for OH, 0.9355
    uE=3.33564*1e-30;
    dE = 1.67;
    h=6.62607*1e-34;
    hbar = h/2/pi;
    c=2.99792458*1e8;
    LD = h*(1.667358e9);


    Bq = [(Bx-1i*By)/sqrt(2), Bz, -(Bx+1i*By)/sqrt(2)];
    Eq = [(Ex-1i*Ey)/sqrt(2), Ez, -(Ex+1i*Ey)/sqrt(2)];

    He = zeros(8);
    Hb = zeros(8);
    H0 = diag(LD*[-ones(1,4) , ones(1,4)]/2);

    for m = -3/2:1:3/2
    for mm = -3/2:1:3/2
    for e = -1:2:1
    for ee = -1:2:1
        te = 0; tb = 0;
        for q = -1:1:1
            te = te + coupling(3/2,m,3/2,e,q,3/2,mm,3/2,ee,-1)*Eq(-q+2);
            tb = tb + coupling(3/2,m,3/2,e,q,3/2,mm,3/2,ee,1)*Bq(-q+2);
        end
        He(m+5/2+2*e+2,mm+5/2+2*ee+2) = te*uE*dE;
        Hb(m+5/2+2*e+2,mm+5/2+2*ee+2) = tb*uB*2*g/.8;

    end
    end
    end
    end

    H = H0 + Hb + He;

end

    function C = coupling(j,m,o,e,q,jj,mm,oo,ee,be)
        term1 = (1+be*e*ee*(-1)^(j+jj+2*abs(o)))/2;
        term2 = sqrt(2*j+1)*sqrt(2*jj+1);
        term3 = Wigner3j([j,1,jj],[-m,q,mm]);
        term4 = Wigner3j([j,1,jj],[-o,0,oo]);
        term5 = (-1)^(m-o);
        term6 = (-1)^q;
        C = term1*term2*term3*term4*term5*term6;
    end

end