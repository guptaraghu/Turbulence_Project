clear all; clc;
%% Data
Hot_wire=dlmread('hotwire.dat');
t=Hot_wire(:,1); l=length(t);
t0=0:2.5e-4:t(end); t0=t0';
dt=t(2)-t(1);
u=Hot_wire(:,2);

%% Calculation
d=trapz(t,u); 
U=d/(t(end)-t(1)); %Avg velocity
U_avg=zeros(l,1)+U; %Avg vel vector
S=u-U_avg; %fluctuating velocity
U_rms=sqrt((1/l)*sum(S.^2)); %Urms


%% Correlation function Calculation
for i=1:1000
    for j=1:l-i+1
          B(i,j)=S(j)*S(i+j-1); 
    end
end
for i=1:1000
   R(i)=(sum(B(i,:))/(l-i+1))/U_rms^2; R=R';
   Y=R(i)*t0(i);
end
Y=sum(Y); %Integral Time-scale
fprintf('Integral time scale= %ds\n',Y)

%Correlation fn plot 
plot(t0(1:20),R(1:20));
xlabel('interval dt');
ylabel('Normalised Correlation function R');
title('Correlation Plot');
grid on
hold on

%% TAYLOR MICRO-SCALE
point1x=0; point1y=1;
point2x=2*dt; point2y=R(3);
point3x=-2*dt; point3y=R(3);

A=[(point1x)^2 point1x 1;
    (point2x)^2 point2x 1;
    (point3x)^2 point3x 1];

B=[point1y;
    point2y;
    point3y];
X=A\B;
a=X(1) ; b= X(2) ; c= X(3);
parabola_coeffs=[a b c];
r=roots(parabola_coeffs);
lambda_val=max(r);  %lambda is the taylor microscale.
x_parabola=linspace(0,lambda_val,100);
y_parabola=a*x_parabola.^2+b*x_parabola+c;

%plot
plot(x_parabola,y_parabola);      
ylim([0 1])   
fprintf('Taylor micro-scale= %ds\n',lambda_val);
 






