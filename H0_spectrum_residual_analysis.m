ind = train_inds(1);
coord = 2;

s_dt = subsamp(1);
T_fin = T_final(1)/length(t);
tobs = t(1:s_dt:end*min(1,T_fin));
dt = mean(diff(tobs));

w_fd = 3; toggle_noise=1;
try
    zdat = [Uobj.Uobs{:}] - (1-toggle_noise)*[Uobj.noise{:,1}];
catch
    zdat = [Uobj.Uobs{:}];
end
xH = arrayfunvec(zdat',rhs_H_true,1)';
c = fdcoeffF(1,0,(-w_fd:w_fd)*dt);
zdot = conv(zdat(:,coord),flipud(c(:,end)),'valid');
res = zdot-xH(w_fd+1:end-w_fd,coord);

m_tf = 5*ppTf;%tf{1}{coord}.rads;%
phi = phifun;%optTFcos(2,0);optTFpoly(1,0);
c = phi((-m_tf:m_tf)/m_tf);
cp = dphi(phi,1); 
cp = [0 cp((-m_tf+1:m_tf-1)/m_tf)/(m_tf*dt) 0];
cp = cp/norm(c,1);
c = c/norm(c,1);
zdot_c = conv(zdat(:,coord),cp,'valid');
xH_c = conv(xH(:,coord),c,'valid');
% zdot_c = conv(zdot,c,'valid');
% xH_c = conv(xH(w_fd+1:end-w_fd,coord),c,'valid');

clf;
subplot(3,2,1:2)
plot(tobs,zdat(:,coord),tobs(1:2*m_tf+1),c/max(c),'linewidth',2)
title('data vs test function')

subplot(3,2,3)
plot(zdot,'ro','DisplayName','z^\prime','linewidth',3); 
hold on;
plot(xH(w_fd+1:end-w_fd,coord),'k-','DisplayName','X_{H_0}(z)','linewidth',3); 
hold off;
title('finite difference DZ vs X_{H_0}(Z)')
fprintf('agreement with H_0: %1.3f\n',norm(res)/norm(zdot))
legend('location','best')
yl = get(gca,'ylim');
set(gca,'ylim',yl,'fontsize',14)
grid on

subplot(3,2,4)
plot(zdot_c,'r','DisplayName','\phi^\prime*z','linewidth',3); 
hold on;
plot(xH_c,'k--','DisplayName','\phi*X_{H_0}(z)','linewidth',3); 
hold off;
title('\phi*X_{H_0}(Z) vs \phi^\prime*Z')
fprintf('agreement with H_0 under phi: %1.3f \n',norm(zdot_c-xH_c)/norm(zdot_c))
legend('location','best')
set(gca,'ylim',yl,'fontsize',14)
grid on

subplot(3,2,5)
semilogy(abs(fft(zdot)),'displayname','FFT(z^\prime)','linewidth',2)
hold on
semilogy(abs(fft(xH(:,coord))),'displayname','FFT(X_{H_0}(z))','linewidth',2)
title('Fourier transform comparison')
hold off
xlim([1 length(zdot)/2])
legend('location','best')
yl = get(gca,'ylim');
set(gca,'ylim',yl,'fontsize',14)
grid on

subplot(3,2,6)
semilogy(abs(fft(xH_c)),'displayname','FFT(\phi*X_{H_0}(z))','linewidth',2)
hold on;
semilogy(abs(fft(zdot_c)),'displayname','FFT(\phi^\prime*z)','linewidth',2)
title('Fourier transform comparison')
xlim([1 length(zdot_c)/2])
hold off;
legend('location','best')
set(gca,'ylim',yl,'fontsize',14)
grid on