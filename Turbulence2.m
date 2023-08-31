clear all; clc;
%% Data
PIV=dlmread('PIVdata.dat');
x=PIV(:,1);
y=PIV(:,2);
u=PIV(:,3);
v=PIV(:,4);
l=length(u);
L=l/400;

%% Ensemble average
en_av=zeros(L,2);
for j=1:1:L
   u0 =0;
    v0=0;
    for i=j:L:l
        u0=u0+u(i);
        v0=v0+v(i);

        en_av(j,1)=u0;
        en_av(j,2)=v0;
    
    end
end
en_av= en_av/400;

fluc = zeros(l,2);
for i=0:1:399
    for j=1:1:14400
        fluc(j+i*14400,1) = u(j+i*14400) - en_av(j,1);
        fluc(j+i*14400,2) = v(j+i*14400) - en_av(j,2);
    end
end

% for x=15.3mm
en_av_x15 = zeros(48000,2);
for i=50:120:l
    m = ((i-50)/120)+1;
    en_av_x15(m,1) = fluc(i,1);
    en_av_x15(m,2) = fluc(i,2);

end

%Reynolds stress Calculation
[uu,uv,vu,vv] = deal(zeros(120,1));
for j=1:1:120
    uu0=0;uv0=0;vu0=0;vv0=0;
    for i=0:120:47880
    uu0 = uu0 + en_av_x15(i+j,1)*en_av_x15(i+j,1);
    uv0 = uv0 + en_av_x15(i+j,1)*en_av_x15(i+j,2);
    vu0 = vu0 + en_av_x15(i+j,2)*en_av_x15(i+j,1);
    vv0 = vv0 + en_av_x15(i+j,2)*en_av_x15(i+j,2);
    end
    uu(j,1) = -1*uu0/400;
    uv(j,1) = -1*uv0/400;
    vu(j,1) = -1*vu0/400;
    vv(j,1) = -1*vv0/400;
end
for i=1:120:14400
    Y(((i-1)/120)+1) = y(14401-i);
end

%%Plot
figure(1)
plot(Y,uu,'LineWidth',2,'DisplayName','Reynolds stress uu')
xlabel('Y');
ylabel('stress');
grid on
hold on
plot(Y,uv,'LineWidth',2,'DisplayName','Reynolds stress uv')
plot(Y,vv,'LineWidth',2,'DisplayName','Reynolds stress vv')
%legend('stress uu','stress uv','stress vv')
legend('Location','southwest','FontSize',12);
hold off

%% Correlation function
y_vector = unique(sort(PIV(:,2)));
x_vector = unique(sort(PIV(:,1)));
snaps = 400; 
xpoint=15.3;
data_points = size(PIV,1)/snaps;
x_req_indicies = find(abs((PIV([1:data_points],1)-xpoint))<0.01);
inputData_cell = {};

for index = 1:snaps
    inputData_cell(index,1) = mat2cell(PIV(data_points*(index-1)+1:(data_points*index),:),[data_points]);
end

y_point = 3.9;
y_indicies= find(abs((PIV([1:data_points],2)-y_point))<=0.01);
y_index = intersect(x_req_indicies,y_indicies);

sum_u_xy = 0;
for snap_index = 1:snaps
    sum_u_xy = sum_u_xy + inputData_cell{snap_index,1}(y_index,3);
end
u_mean_xy = sum_u_xy/snaps;

u_rss_sum = 0;
for snap_index = 1:snaps
    u_rss_sum = u_rss_sum + ((inputData_cell{snap_index,1}(y_index,3)-u_mean_xy)^2);
end

u_rms_xy = sqrt(u_rss_sum/snaps);
Ruu = zeros(data_points,1);
u_velocity = zeros(data_points,snaps);
u_mean_velo = zeros(data_points,1);

for n_data_point = 1:data_points
    
    velo_sum = 0;
    
    for snap_index = 1:snaps
        velo_sum = velo_sum + inputData_cell{snap_index,1}(n_data_point,3);
    end
    
    u_mean_velo(n_data_point) = velo_sum/snaps;
    
end

for n_data_point = 1:data_points
    
    for snap_index = 1:snaps
        u_velocity(n_data_point,snap_index) = inputData_cell{snap_index,1}(n_data_point,3) - u_mean_velo(n_data_point,1);
    end

end

for n_data_point = 1:data_points
    
    Ruu_sum = 0;
    for snap_index = 1:snaps
        Ruu_sum = Ruu_sum + u_velocity(n_data_point,snap_index)*u_velocity(y_index,snap_index);
    end
    
    Ruu(n_data_point) = Ruu_sum/(snaps*(u_rms_xy^2));
end

Ruu_corr_data = zeros(round(sqrt(data_points)),round(sqrt(data_points)));
Ruu_corr_data = transpose(reshape(Ruu,size(Ruu_corr_data,1),size(Ruu_corr_data,1)));

%% Plotting the Ruu correlation contour
figure('Color','white','Position',[250 70 950 700]);
contour(x_vector,y_vector,flipud(Ruu_corr_data),400);
set(gca(),'FontSize',12);
xlabel('X \rightarrow');
ylabel('Y \rightarrow');
title('Correlation contour at x_{0} \approx 15.3 mm and y_{0} \approx 3.9 mm');
colorbar();

%% Calculation of Taylor longitudinal microscale

x_max_indicies = find(PIV(1:data_points) == max(x_vector));

for x = 1:length(x_max_indicies)
    
    if x_max_indicies(x) > y_index
        x_max_req_index = x_max_indicies(x);
        break;
    else
        continue;
    end
end

%Calculate Longitudinal correlation function
corr_long_data = zeros(length(y_index:x_max_req_index),1);
del_x_data = zeros(length(y_index:x_max_req_index),1);
counter = 1;
for index = y_index:x_max_req_index
    del_x_data(counter,1) = inputData_cell{1,1}(index,1) - inputData_cell{1,1}(y_index,1);
    Rxx_sum = 0;
    
    for snap_index = 1:snaps
        Rxx_sum = Rxx_sum + u_velocity(index,snap_index) * u_velocity(y_index,snap_index);
    end
    corr_long_data(counter,1) = Rxx_sum/(snaps*u_rms_xy^2);
    counter = counter+1;
end

%Defining data points for parabola fitting
point1_x = 0;
point1_y = 1;

point2_x = del_x_data(2);
point2_y = corr_long_data(2);

point3_x = -del_x_data(2);
point3_y = corr_long_data(2);

%Parabola is given in matrix form as [A]*[X].^2 + [B].*[X] + C = Y

%Defining [A] matrix
A_long=[(point1_x)^2 point1_x 1;
    (point2_x)^2 point2_x 1;
    (point3_x)^2 point3_x 1];

%Defining [B] matrix
B_long=[point1_y;
    point2_y;
    point3_y];

%Calculating the coefficients
X_long=A_long\B_long;

%Assigning the coefficients
a=X_long(1) ; b= X_long(2) ;c=X_long(3);
parabola_coeffs_long=[a b c];

%The roots of the parabola will give where the parabola will intersect with
%x axis which will give us the Taylor microscale
x_crossing_val=roots(parabola_coeffs_long);

%Taylor microscale
lambda_long_val=max(x_crossing_val);  %lambda is the taylor microscale.
x_parabola_long=linspace(-lambda_long_val,lambda_long_val,100);
y_parabola_long=a*x_parabola_long.^2+b*x_parabola_long+c;

%Plotting longitudinal correlation function and microscale
figure('Color','white','Position',[250 70 950 700]);
p_handle = plot(vertcat(-flipud(del_x_data),del_x_data),...
                vertcat(flipud(corr_long_data),(corr_long_data)));
set(p_handle,'LineWidth',2,'DisplayName','Longitudinal Correlation function');
hold on;
set(gca(),'FontSize',12);
plot(x_parabola_long,y_parabola_long,'LineWidth',2,'DisplayName','Fitted - Parabola curve');
xlabel('X \rightarrow');
ylabel('Longitudinal correlation function)');
title('Longitudinal correlation function & Longitudinal Taylor microscale estimation');
grid on;
legend();

%% Calculation for Taylor transverse microscale
del_y_data = [];
corr_trans_data = [];

%Calculate transverse correlation function
for index_xy = 1:length(x_req_indicies)
    
    if x_req_indicies(index_xy) >= y_index
       break; 
    end
    
    del_y_data(index_xy,1) = inputData_cell{1,1}(x_req_indicies(index_xy),2)-inputData_cell{1,1}(y_index,2);
    
    Rxy_sum = 0;
    for snap_index = 1:snaps
        Rxy_sum = Rxy_sum + u_velocity(y_index,snap_index)*u_velocity(x_req_indicies(index_xy),snap_index);
    end
   
    corr_trans_data(index_xy,1) = Rxy_sum/(snaps*u_rms_xy^2);
end
del_y_data = flipud(del_y_data);
corr_trans_data = flipud(corr_trans_data);

%Defining data points for parabola fitting
point1_x = 0;
point1_y = 1;

point2_x = del_y_data(2);
point2_y = corr_trans_data(2);

point3_x = -del_y_data(2);
point3_y = corr_trans_data(2);

%Parabola is given in matrix form as [A]*[X].^2 + [B].*[X] + C = Y

%Defining [A] matrix
A_trans=[(point1_x)^2 point1_x 1;
    (point2_x)^2 point2_x 1;
    (point3_x)^2 point3_x 1];

%Defining [B] matrix
B_trans=[point1_y;
    point2_y;
    point3_y];

%Calculating the coefficients
X_trans=A_trans\B_trans;

%Assigning the coefficients
a=X_trans(1) ; b= X_trans(2) ;c=X_trans(3);
parabola_coeffs_long=[a b c];

%The roots of the parabola will give where the parabola will intersect with
%x axis which will give us the Taylor microscale
x_crossing_val=roots(parabola_coeffs_long);

%Taylor microscale
lambda_trans_val=max(x_crossing_val);  %lambda is the taylor microscale.
x_parabola_trans=linspace(-lambda_trans_val,lambda_trans_val,100);
y_parabola_trans=a*x_parabola_trans.^2+b*x_parabola_trans+c;

figure('Color','white','Position',[250 70 950 700]);
p_handle = plot(vertcat(-flipud(del_y_data),del_y_data),...
                vertcat(flipud(corr_trans_data),(corr_trans_data)));
set(p_handle,'LineWidth',2,'DisplayName','Transverse Correlation function');
hold on;
set(gca(),'FontSize',12);
plot(x_parabola_trans,y_parabola_trans,'LineWidth',2,'DisplayName','Fitted - Parabola curve');
xlabel('Y \rightarrow');
ylabel('Transverse correlation function)');
grid on
title('Transverse correlation function & Transverse Taylor microscale estimation');






