% For OH:
A = -4168639.13;
B = 555660.97;
g = -3574.88;
p = 7053.09846;
q = -1159.99165;

% For OD:
A = -4165275.23
B = 296312.0
g = -10000; % I don't know what this one is. But it doesn't change the lambda doublet TOO much.
p = 3765.112;
q = -329.38;

% Conclusion: Lambda Doublet for OD ~ 300 MHz.

m11 = A/2 + 2*B;
m22 = -A/2 + 4*B - g + p + 2*q;
m12 = -sqrt(3)*(B-g/2+q);

m33 = A/2 + 2*B;
m44 = -A/2 + 4*B - g - p - 2*q;
m34 = -sqrt(3)*(B-g/2-q);



eig([m11 m12 ; m12 m22])-eig([m33 m34 ; m34 m44])