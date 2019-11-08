load('DMreqn0exposure0.1and1.mat')
save1 = sig5SiSave;
load('DMreqn2exposure0.1first99correct.mat')
save2 = sig5Si;
load('DMreqn2exposure0.1restcorrect.mat')
save3 = sig5Si;
save2(100:109) = save3(1:10);
load('DMreqn2exposure1first79correct.mat')
save4 = sig5Si;
load('DMreqn2exposure1restcorrect.mat');
save5 = sig5Si;
save4(80:109) = save5(1:30);

save1(3,:) = save2;
save1(4,:) = save4;

min(save1(1,:))
min(save1(2,:))
min(save1(3,:))
min(save1(4,:))

save1 = log10(save1)
plot(1:109,save1(1,:))
hold on
plot(1:109,save1(2,:))
plot(1:109,save1(3,:))
plot(1:109,save1(4,:))



