close all
v0Array = linspace(200,260,51);                                       
vEArray = linspace(200,280,51);                     
vescArray = linspace(500,700,51); 

load('FDM_1_vesc')
save1 = sigSi_5;
load('FDM_1_vE')
save2 = sigSi_5;
load('FDM_1_v0')
save3 = sigSi_5;
load('FDM_q2_v0.mat')
save4 = sigSi_5;
load('FDM_q2_vE.mat')
save5 = sigSi_5;
load('FDM_q2_vesc.mat')
save6 = sigSi_5;

figure(1)
plot(v0Array,save3)
hold on
plot(v0Array,save4)
axis([200 260 0 9*10^-40])

figure(2)
plot(vEArray,save2)
hold on
plot(vEArray,save5)
axis([200 280 0 9*10^-40])

figure(3)
plot(vescArray,save1)
hold on
plot(vescArray,save6)
axis([500 700 0 9*10^-40])