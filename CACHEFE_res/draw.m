close all;clear all;
cd ~/Documents/research/cacti65/code/APA' zqgao'/try3/CACHEFE_res/;
cells=[1:8,10:3:16,32:16:64];

Fp1=load('./APAdemo_safail__res.txt');nCend=8;
figure(1);
for j=1:8
    plot(cells(cells>=j-1),Fp1(j,cells>=j-1),'LineWidth',2);
    hold on;
end
title('SA Failure');

Fp2=load('./APAdemo_writefail__res.txt');nCend=8;
figure(2);
for j=1:8
    plot(cells(cells>=j-1),Fp2(j,cells>=j-1),'LineWidth',2);
    hold on;
end
title('Write Failure');


Fp3=load('./APAdemo_readfail__res.txt');nCend=8;
figure(3);
for j=1:8
    plot(cells(cells>=j-1),Fp3(j,cells>=j-1),'LineWidth',2);
    hold on;
end
title('Read Failure');
