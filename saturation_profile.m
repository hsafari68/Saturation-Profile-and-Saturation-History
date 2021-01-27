%% Start
clc;
clear;
close all;

%% Import necessary data
x = xlsread('Position.xlsx'); % Import position matrix
t = xlsread('Time.xlsx'); % Import time matrix in seconds
xray = xlsread('Xray.xlsx'); % Import low-high sequences of xray hits
GL = xlsread('Gas_L.xlsx'); % Import Gas low energy xray
GH = xlsread('Gas_H.xlsx'); % Import Gas high energy xray
OL = xlsread('Oil_L.xlsx'); % Import Oil low energy xray
OH = xlsread('Oil_H.xlsx'); % Import Oil high energy xray
WL = xlsread('Water_L.xlsx'); % Import Water low energy xray
WH = xlsread('Water_H.xlsx'); % Import Water high energy xray

%% Interpolation
% Six cubic interpolations for caliberation data
pgl = spline (GL(:,1),GL(:,2)); % Low energy gas xray 
pgh = spline (GH(:,1),GH(:,2)); % High energy gas xray
pol = spline (OL(:,1),OL(:,2)); % Low energy oil xray
poh = spline (OH(:,1),OH(:,2)); % High energy oil xray
pwl = spline (WL(:,1),WL(:,2)); % Low energy water xray
pwh = spline (WH(:,1),WH(:,2)); % High energy water xray


% Two bicubic interpolations for low-high sequences of X-ray hits
% Low and high X-ray time and hits should be separated first as well as
% their positions

tL = zeros(size(t,1),16);
tH = zeros(size(t,1),16);
L = zeros(size(t,1),16);
H = zeros(size(t,1),16);
xL = zeros(size(t,1),16);
xH = zeros(size(t,1),16);

for i =1:16
    xL(:,i) = x(:,2*i-1); 
    xH(:,i) = x(:,2*i);
    tL(:,i) = t(:,2*i-1); 
    tH(:,i) = t(:,2*i);
    L(:,i) = xray(:,2*i-1);
    H(:,i) = xray(:,2*i);
    
end

% In order to obtain saturation profiles, the low and high X-ray hits
% need to be interpolated at positions where the caliberation data are
% available.


% LGinterp = interp2 (x,tL,L,GL(:,1),tL);
% HGinterp = interp2 (x,tH,H,HL(:,1),tH);
% LOinterp = interp2 (x,tL,L,OL(:,1),tL);
% HOinterp = interp2 (x,tH,H,OH(:,1),tH);
% LWinterp = interp2 (x,tL,L,WL(:,1),tL);
% HWinterp = interp2 (x,tH,H,WH(:,1),tH);

% Interpolating the caliberation data at the points where low and high
% X-ray hits have been recorded

LGinterp = ppval (pgl,x);
HGinterp = ppval (pgh,x);
LOinterp = ppval (pol,x);
HOinterp = ppval (poh,x);
LWinterp = ppval (pwl,x);
HWinterp = ppval (pwh,x);

% Interpolating the low and high X-ray to have their values at all times
% and spaces

Htotal = griddata (xH,tH,H,x,t);
Ltotal = griddata (xL,tL,L,x,t);


%% Solving the linear system of equations to obtain saturation profile
% So*logHo+Sg*logHg+Sw*logHw = Log H*
% So*logLo+Sg*logLg+Sw*logLw = Log L*
% So+Sg+Sw = 1

% Please note that this program calculates the saturation quantities at
% each point by solving the three linear equations separately for each
% position and time

So = zeros(size(t));
Sg = zeros(size(t));
Sw = zeros(size(t));

for ii = 1: size(t,1)
    for jj = 1:size(t,2)
        
        A = [log(LOinterp(ii,jj)) log(LGinterp(ii,jj)) log(LWinterp(ii,jj))
             log(HOinterp(ii,jj)) log(HGinterp(ii,jj)) log(HWinterp(ii,jj))
             1 1 1
             ];
         B = [log(Ltotal(ii,jj)) log(Htotal(ii,jj)) 1]';
         S = linsolve(A,B);
         So (ii,jj) = S(1);
         Sg (ii,jj) = S(2);
         Sw (ii,jj) = S(3);
        
 
    end
end

%% Visualization

subplot(2,1,1)
mesh(xL,tL,L);
xlabel('\bfPosition, mm');
ylabel('\bfTime, s');
zlabel('\bfX-ray, hits/s');
title('\bfLow energy X-ray');
subplot(2,1,2)
mesh(xH,tH,H);
xlabel('\bfPosition, mm');
ylabel('\bfTime, s');
zlabel('\bfX-ray, hits/s');
title('\bfHigh energy X-ray');

figure(2)
subplot(2,2,1)
plot(x,LOinterp,'-b','LineWidth',1.5);
xlabel('\bfPosition, mm');
ylabel('\bfLow energy X-ray, hits/s');
subplot(2,2,[2,4])
mesh(x,t,So);
xlabel('\bfPosition, mm');
ylabel('\bfTime, s');
zlabel('\bfOil saturation (S_{o})');
subplot(2,2,3)
plot(x,HOinterp,'-b','LineWidth',1.5);
xlabel('\bfPosition, mm');
ylabel('\bfHigh energy X-ray, hits/s');

figure(3)
subplot(2,2,1)
plot(x,LGinterp,'-b','LineWidth',1.5);
xlabel('\bfPosition, mm');
ylabel('\bfLow energy X-ray, hits/s');
subplot(2,2,[2,4])
mesh(x,t,Sg);
xlabel('\bfPosition, mm');
ylabel('\bfTime, s');
zlabel('\bfGas saturation (S_{g})');
subplot(2,2,3)
plot(x,HGinterp,'-b','LineWidth',1.5);
xlabel('\bfPosition, mm');
ylabel('\bfHigh energy X-ray, hits/s');


figure(4)
subplot(2,2,1)
plot(x,LWinterp,'-b','LineWidth',1.5);
xlabel('\bfPosition, mm');
ylabel('\bfLow energy X-ray, hits/s');
subplot(2,2,[2,4])
mesh(x,t,Sw);
xlabel('\bfPosition, mm');
ylabel('\bfTime, s');
zlabel('\bfWater saturation (S_{w})');
subplot(2,2,3)
plot(x,HWinterp,'-b','LineWidth',1.5);
xlabel('\bfPosition, mm');
ylabel('\bfHigh energy X-ray, hits/s');

