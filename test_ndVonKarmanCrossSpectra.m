close all
L0 = 340;
H = 150;
x_ad = linspace(0,1,1001);

t_c    = 3.87 /12;
fmin = log10(1e-3);
fmax = log10(10);
nf   = 501;
f_i  = logspace(fmin,fmax,nf)';
kapa = t_c * f_i;
[Su,coh] = ndVonKarmanCrossSpectra(kapa,x_ad);

Su_kaimal   = kaimalSpectra(kapa,1,1,1);

[xMesh,fMesh] = meshgrid(x_ad,f_i);
figure(1)
contour(xMesh,fMesh,coh,'ShowText','on')
xlabel('$x/L_{0}$ [-]')
ylabel('f [Hz]')

% figure(2)
% plot(x_ad,coh(1,:),'-b');hold on
% plot(x_ad,coh(100,:),'-k');hold on
% plot(x_ad,coh(200,:),'-r');hold on
% plot(x_ad,coh(300,:),'-y');hold on
% plot(x_ad,coh(400,:),'-m');hold on
% plot(x_ad,coh(500,:),'-c');hold on
% plot(x_ad,linspace(1,1,1000),'--k');hold on
% xlabel('$x/L_{0}$ [-]')
% ylabel('Coh [-]')
% legend('f = 0.001 Hz','f = 0.0062 Hz','f = 0.0391 Hz','f = 0.2466 Hz','f = 1.556 Hz','f = 9.8175 Hz')

figure(2)
plot(f_i,coh(:,1),'k-'); hold on; 
plot(f_i,coh(:,201),'b-'); hold on; 
plot(f_i,coh(:,401),'r-'); hold on; 
plot(f_i,coh(:,601),'y-'); hold on; 
plot(f_i,coh(:,801),'g-'); hold on; 
plot(f_i,coh(:,1001),'c'); hold on; 
xlabel('f [Hz]')
ylabel('Coh [-]')

figure(3)
plot(x_ad,coh(1,1),'k-'); hold on; 
plot(x_ad,coh(101,:),'b-'); hold on; 
plot(x_ad,coh(201,:),'r-'); hold on; 
plot(x_ad,coh(301,:),'y-'); hold on; 
plot(x_ad,coh(401,:),'g-'); hold on; 
plot(x_ad,coh(501,:),'c'); hold on; 
xlabel('x/L [-]')
ylabel('Coh [-]')

figure(4)
loglog(f_i,kapa.*Su_kaimal,'k--'); hold on; 
loglog(f_i,kapa.*Su(:,1),'k-'); hold on; 
loglog(f_i,kapa.*Su(:,201),'b-'); hold on; 
loglog(f_i,kapa.*Su(:,401),'r-'); hold on; 
loglog(f_i,kapa.*Su(:,601),'y-'); hold on; 
loglog(f_i,kapa.*Su(:,801),'g-'); hold on; 
loglog(f_i,kapa.*Su(:,1001),'c'); hold on; 
xlabel('f [Hz]')
ylabel('f Su [Hz]')