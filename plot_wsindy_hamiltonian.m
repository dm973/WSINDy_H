warning('off','MATLAB:legend:IgnoringExtraEntries')

subplot(3,4,1:2)
plot([b_plot G_plot*w_plot G_plot*w_true])
legend({'b','Gw','Gw_true'},'box','off')
title(['||Gw-b||/||b||=',num2str(norm(b_plot-G_plot*w_plot)/norm(b_plot)),'; ||w-w^*||/||w^*||=',num2str(norm(w_true-w_plot)/norm(w_true))])

subplot(3,4,3:4)
plot(1:J,w_plot,'x',1:J,w_true,'o','linewidth',2,'markersize',10)
grid on
legend({'w','w^*'},'location','best')
title(['TPR=',num2str(tpscore(w_plot,w_true)),'; sparsity=',num2str(length(find(w_plot))),'/',num2str(size(G_plot,2))])
xlabel('Library index')

subplot(3,4,5:6)
x0_foo = num2cell(x0_reduced);
C = Hfun(x0_foo{:});
x0_foo = num2cell(x0_true);
C_true = Hfun_true(x0_foo{:});

if ~isequal(plot_style,'Henon-Heiles')
    bndsQ = [min(arrayfunvec(xobs_plot(:,1:2:end)',@min,2)) max(arrayfunvec(xobs_plot(:,1:2:end)',@max,2))];
    bndsP = [min(arrayfunvec(xobs_plot(:,2:2:end)',@min,2)) max(arrayfunvec(xobs_plot(:,2:2:end)',@max,2))];
    Np = 100;   Nq = 100;
%     wq = range(bnds(1,:))*range_fac; wp = range(bnds(2,:))*range_fac;
    wq = range(bndsQ)*range_fac;    wp = range(bndsP)*range_fac;
    [qq,pp]=meshgrid(linspace(-wq+bndsQ(1),bndsQ(2)+wq,Nq),linspace(-wp+bndsP(1),bndsP(2)+wp,Np));
    if nstates == 2
        F = Hfun(qq,pp);
    elseif nstates == 4
        F = Hfun(qq,pp,qq,pp);
    else
        F = 0*qq;
    end
    cvs = linspace(min(F(:)),max(F(:)),50);
    [~,h1]=contourf(qq,pp,F,cvs,'EdgeAlpha',0.25);
    hold on
    for i=1:2:nstates
        h2=plot(xobs_plot(:,i),xobs_plot(:,i+1),'w-','linewidth',2);
        hold on
        h3=plot(xH0_learned(:,i),xH0_learned(:,i+1),'k--','linewidth',2);
    end
    hold off
    legend([h1;h2;h3],{'H','obs','sim'})
    xlabel('q'); ylabel('p')
    colorbar
    title(['H_{learned}'])
    subplot(3,4,9:10)
    if nstates == 2
        F_true = Hfun_true(qq,pp)-C_true+C;
    elseif nstates == 4
        F_true = Hfun_true(qq,pp,qq,pp)-C_true+C;
    else
        F = 0*qq;
    end
    contourf(qq,pp,F_true,cvs,'EdgeAlpha',0.25)
    hold on
    for i=1:2:nstates
        h1=plot(xobs_plot(:,i),xobs_plot(:,i+1),'w','linewidth',2);
        hold on
        h2=plot(xH0_true(:,i),xH0_true(:,i+1),'k--','linewidth',2);
    end
    xlabel('q'); ylabel('p')
    colorbar
    legend({'H','obs','sim H_{true}'},'Box','off')
    hold off
    title(['H_{true}. ||H_L-H_T||_2/||H_T||_2=',num2str(norm(F(:)-F_true(:))/norm(F_true(:)))])    


else
    bnds = arrayfunvec(xobs_plot(:,[1 3])',@minmax,2);
    Np = 100; Nq = 100;
    wq = range(bnds(1,:))*range_fac; wp = range(bnds(2,:))*range_fac;
    [qq,pp]=meshgrid(linspace(-wq+bnds(1,1),bnds(1,2)+wq,Nq),linspace(-wp+bnds(2,1),bnds(2,2)+wp,Np));
    F = Hfun(qq,0,pp,0);
    cvs = linspace(min(F(:)),max(F(:)),50);
    [~,h1]=contourf(qq,pp,F,cvs,'EdgeAlpha',0.25);
    hold on
    h2=plot(xobs_plot(:,1),xobs_plot(:,3),'w-','linewidth',2);
    h3=plot(xH0_learned(:,1),xH0_learned(:,3),'k--','linewidth',2);
    hold off
    legend([h1;h2;h3],{'H','obs','sim'})
    xlabel('q'); ylabel('p')
    colorbar
    title(['H_{learned}'])

    subplot(3,4,9:10)
    try
        F_true = Hfun_true(qq,0,pp,0)+C-C_true;
        [~,h1]=contourf(qq,pp,F_true,cvs,'EdgeAlpha',0.25);
        hold on
        h2=plot(xobs_plot(:,1),xobs_plot(:,3),'w-','linewidth',2);
        h3=plot(xH0_true(:,1),xH0_true(:,3),'k--','linewidth',2);
        hold off
        legend([h1;h2;h3],{'H','obs','sim'})
        xlabel('q'); ylabel('p')
        colorbar
        title(['H_{true}. ||H_L-H_T||_2/||H_T||_2=',num2str(norm(F(:)-F_true(:))/norm(F_true(:)))])
    catch
        for i=1:2:nstates
            h1=plot(xobs_plot(:,i),xobs_plot(:,i+1),'-','linewidth',2);
            hold on
            h2=plot(xH0_true(:,i),xH0_true(:,i+1),'--','linewidth',2);
        end
        hold off
        legend([h1;h2],{'obs','sim'})
    end
end

% subplot(3,4,9:10)
% try
%     plot(tplot,Hfun(xobs_plot(:,1),xobs_plot(:,2)),'ro',tplot,Hfun(xtrue_plot(:,1),xtrue_plot(:,2)),'k-')
%     title('H_{learned}(q(t),p(t))')
% catch
%     xlabel('t')
%     h1=plot(tplot,xobs_plot,'-',tplot,xH0_learned,'--','linewidth',2);
%     legend(h1([1 end]),{'obs','sim'})
%     title(num2str(norm(xtrue_plot(:)-xH0_learned(:))/norm(xtrue_plot(:))))
% end

% subplot(3,4,10)
% try
%     plot(tplot,Hfun_true(xobs_plot(:,1),xobs_plot(:,2)),'ro',tplot,Hfun_true(xtrue_plot(:,1),xtrue_plot(:,2)),'k-')
%     title('H_{true}(q(t),p(t))')
% catch
%     xlabel('t')
%     h1=plot(tplot,xobs_plot,'-',tplot,xH0_true,'--','linewidth',2);
%     legend(h1([1 end]),{'obs','true'})
%     title(num2str(norm(xtrue_plot(:)-xH0_true(:))/norm(xtrue_plot(:))))
% end

% xfft = abs(fftshift(fft(xtrue_plot./vecnorm(xH0_learned,1,1))));
% xfft_learned = abs(fftshift(fft(xH0_learned./vecnorm(xH0_learned,1,1))));
xfft = abs(fftshift(fft(xtrue_plot)));
xfft_learned = abs(fftshift(fft(xH0_learned)));
ks = [0:length(tplot)-1]-ceil(length(tplot)/2);
ks_learned = [0:size(xH0_learned,1)-1]-ceil(size(xH0_learned,1)/2);

subplot(3,4,7:8)
semilogy(ks,xfft,'g-','linewidth',2);
hold on
semilogy(ks_learned,xfft_learned,'b--','linewidth',2);
if ~isempty(vfft)
    semilogy(kx,arrayfun(@(i)xfft(kx(i),i),1:nstates),'ko','linewidth',4,'markersize',10);
    semilogy(ks,vfft./max(vfft(:))*max(xfft(:))','r-','linewidth',2);
end
hold off
xlim(ks([ceil(end/2) end]))
% ylim([min(xfft(:))/10 1])
title(['m_t=',num2str(mt)])
lns=get(gca,"Children");
legend(lns(1:nstates:end),{'fft learned x','fft true'},'box','off','location','sw')

try
    xfft_true = abs(fftshift(fft(xH0_true./vecnorm(xH0_true,1,1))));
    subplot(3,4,11:12)
    lns=plot(tplot,xtrue_plot,'g-',tplot(1:size(xH0_true,1)),xH0_true,'k--',tplot(1:size(xH0_learned,1)),xH0_learned,'b-','linewidth',2);
    hold on
    tfs=arrayfun(@(i)plot(tplot(1:2*mt(i)+1),v{i}+i-1,'r-'),1:nstates);
    hold off
%     legend([lns(1:nstates:end);tfs(1)],{'x','x red', 'x learned','\phi'})
    legend([lns(1:nstates:end)],{'x','x red', 'x learned','\phi'})
    try
    title(['true red err=',num2str(norm(xH0_true(:)-xtrue_plot(:))/norm(xtrue_plot(:))),...
        '; learned err=',num2str(norm(xH0_learned(:)-xtrue_plot(:))/norm(xtrue_plot(:))),...
        '; learned-true red=',num2str(norm(xH0_learned(:)-xH0_true(:))/norm(xH0_true(:)))])
    end
catch
end
