clear all;
% File which includes three columns, y U U''
data= load ('velocityprofile1.dat');  % load your velocity profile
%----------------------------------------

% Number of Chebyshev modes
N= 99; N1 = N+1; 
j = (0:N)'; X = cos(pi*j/N);
%-------------------------------------

% Scaling for mapping from -1 to 1
L = data(end,1);
scal=L/2; % Physical domain => [0, L]
Y=X*scal + scal; % Mapping from Physical to Chebyshev domain
%----------------------------------

% Velocity and its second derivative are being evaluated at scaled Y
U = interp1(data(:,1), data(:,2), Y);
U2 = interp1(data(:,1), data(:,3), Y);
fac = 1/scal;
%-------------------------------------

%%For singel value of alpha and Re, for example,
alpha = 0.3;
length(alpha)
Rey = 600;
%%-------------------------------

% For a range of alpha and Reynolds number, for example,
% alpha = 0.1:0.01:1.5;
% length(alpha)
% Rey = 10:0.5:3000;
%-----------------------------------------

% Main program; try to understand before you run it
ci = zeros(length(alpha),length(Rey));
cr = zeros(length(alpha),length(Rey));
for ii = 1:length(alpha)
    ii
    for jj = 1:length(Rey)
        al = alpha(ii);
        R = Rey(jj);
        zi = sqrt(-1); a2 = al^2; a4 = a2^2; er = -200*zi;
        [D0,D1,D2,D3,D4] = Dmat(N); % Read about it in the book.
        D1 = fac*D1;
        D2 = (fac^2)*D2;
        D4 = (fac^4)*D4;
        B = (D2-a2*D0);
        A = (U*ones(1,N1)).*B-(U2*ones(1,N1)).*D0-(D4-2*a2*D2+a4*D0)/(zi*al*R);
        A = [er*D0(1,:); er*D1(1,:); A(3:N-1,:); er*D1(N1,:); er*D0(N1,:)];
        B = [D0(1,:); D1(1,:); B(3:N-1,:); D1(N1,:); D0(N1,:)];
        d = (inv(B)*A);
        [vv, c] = eig(d);                 % eigenvalues are being evaluated using eig function
%         [mxci,I] = max(imag(c));  % useful for contour plots
%         ci(ii,jj)= mxci;
%         cr(ii,jj) =real(c(I)); 

    end;
end
