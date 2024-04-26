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

%%

zdot = Uobj.Uobs{j};
km = 10; N = 1000;
x = fftshift(abs(fft(zdot)));
x = x(1:floor(end/2));
[~,Umax] = max(x);
y = cumsum(x(1:Umax)); y = y/y(end);
z = (1:N)/N;
[a,b,c]=kmeans(arrayfun(@(x) find(y>=x,1),z)',km);
length(x)-sum(arrayfun(@(i)b(i)*length(find(a==i))/length(a),1:km))
plot(arrayfun(@(i)b(i)*length(find(a==i))/length(a),1:km))

%%
j=1;
i=1;
toggle_plot=0;
xobs = permute(Uobj.Uobs{j},[i 1:i-1 i+1:max(2,Uobj.ndims)]);
t = Uobj.grid{i}(:);
T = length(t);
wn = ((0:T-1)-floor(T/2))'*(2*pi)/range(t);
xx = wn(1:ceil(end/2));
NN = length(xx);

Ufft = mean(reshape(abs(fftshift(fft(xobs))),T,[]),2) /sqrt(NN);
Ufft = Ufft(1:ceil(T/2),:);
[~,Umax] = max(Ufft);
cornerdat = cumsum(Ufft(1:Umax));

tstarind = getcorner(cornerdat,xx(1:Umax),tplot);
max(ceil(length(Ufft)-tstarind+1),1)
% Uobj.ks(j,i) = max(ceil((length(Ufft)-tst

%%

y = res;
N = length(y);

m1 = 25;
b1 = convmtx(phifun((-m1:m1)/m1)',N);
b1 = b1(2*m1+1:end-2*m1,:);

m2 = 50;
b2 = convmtx(phifun((-m2:m2)/m2)',N);
b2 = b2(2*m2+1:end-2*m2,:);

G = eye(N) - y(:)*y(:)'/norm(y)^2;
D = [zeros(m1,N);b1*dftmtx(N);zeros(m1,N)];
K = N-2*m2-1;
A = G*D(:,1:K);

W = A \ b2';
E = real(A*W);
% E = E-E(1,:);
norm(E'*y(:))/norm(y(:))

for i=1:N-2*m
    subplot(3,1,1)
    plot(E(:,i),'o-')
    hold on
    plot(y,'k--')
    hold off
    title(dot(E(:,i),y)/norm(y(:)))
    subplot(3,1,2)
    plot(b2(i,:),'o-')
    subplot(3,1,3)
    semilogy(imag(fft(E(:,i))),'o-')
    hold on
%     semilogy(imag(fft(b(:,i))),'o-')
    semilogy(imag(fft(y)),'k-x')
    title(dot(fft(y),fft(E(:,i)))/norm(y(:)))
    hold off
    drawnow
end
