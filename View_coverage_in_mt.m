%% coverage in mt
dr = '~/Dropbox/Boulder/research/data/WSINDy_H_data/WSINDy_H_materials/';
for expnum = 1:4
for ii=5
mts_res = 80;

if expnum==1
load([dr,'expnum',num2str(expnum),'_mtsweep.mat'],'Herr_all','Xerr_all','X0err_all','X0_trueerr_all','tpr_all','mt_test_all','epz','mts','ppTf','train_inds')
plot_inds =[8:2:30];%exp1
epn = [0.01 0.05];
elseif expnum==3
load([dr,'expnum',num2str(expnum),'_mtsweep.mat'],'Herr_all','Xerr_all','X0err_all','X0_trueerr_all','tpr_all','mt_test_all','epz','mts','ppTf','train_inds')
plot_inds = 1:10;%exp3
epn = [0.01 0.03];
elseif expnum==4
load([dr,'expnum',num2str(expnum),'_mtsweep.mat'],'Herr_all','Xerr_all','X0err_all','X0_trueerr_all','tpr_all','mt_test_all','epz','mts','ppTf','train_inds')
plot_inds = [1:2:12];%exp4
epn = [0.01 0.03];
end

dats = {Herr_all(:,:,end-4:end),...
    Xerr_all(:,:,end-4:end),...
    X0err_all(:,:,end-4:end),...
    X0_trueerr_all(:,:,end-4:end),...
    tpr_all(:,:,end-4:end)};

for i=length(epn)
    mt_learned = min(squeeze(mt_test_all(plot_inds,1,epn(i)==epz))/ppTf*2,max(mts));

    dat = reshape(dats{ii},length(train_inds),[],5);
    dat = min(squeeze(dat(plot_inds,:,epz==epn(i))),1);
    dat_learned = diag(dat(:,ceil(mt_learned)));
    disp([ expnum;find(dat_learned<1)])
    figure(ii);clf;

    imagesc(dat);set(gca,'Ydir','normal')
    hold on
    plot(mt_learned,1:length(plot_inds),'k-o','linewidth',4);hold off
    ylim([1 length(plot_inds)])
    xlim([1 mts_res])
    
    if ii==5
        clim([0.5 1])
        colorbar('ticklabelinterpreter','latex')
        colormap(gray)
    else
        clim([0.01 1])
        set(gca,'ColorScale','log')
        colorbar('ticklabelinterpreter','latex','XTick', 10.^(-2.5:0.5:0))
        colormap(cool)
    end
    set(gca,'fontsize',16,'ticklabelinterpreter','latex')
    xlabel('$\sigma_{\phi f}$','interpreter','latex')
    ylabel('$z(0)$ index','interpreter','latex') 
    drawnow
    saveas(gcf,['~/Desktop/expnum',num2str(expnum),'_ep',num2str(epn(i)),'_stat',num2str(ii),'_mt.png'])
end

if expnum==2
    mQs=[1 5];
    for i=length(mQs)
        mQ = mQs(i);
        load([dr,'expnum',num2str(expnum),'_mQ',num2str(mQ),'_mtsweep.mat'],'Herr_all','Xerr_all','X0err_all','X0_trueerr_all','tpr_all','mt_test_all','epz','mts','ppTf','train_inds')
        plot_inds = [7:5:48];%exp2
        ep = 0.05;
        
        dats = {Herr_all(:,:,end-4:end),...
            Xerr_all(:,:,end-4:end),...
            X0err_all(:,:,end-4:end),...
            X0_trueerr_all(:,:,end-4:end),...
            tpr_all(:,:,end-4:end)};
        
        mt_learned = min(squeeze(mt_test_all(plot_inds,1,ep==epz))/ppTf*2,max(mts));
    
        if ii==5
            dat = reshape(dats{ii},length(train_inds),[],5);
            clim([0.75 1])
        else
            dat = reshape(dats{ii},length(train_inds),[],5);
        end
        dat = min(squeeze(dat(plot_inds,:,epz==ep)),1);
        dat_learned = diag(dat(:,ceil(mt_learned)));
        disp([ expnum;find(dat_learned<1)])
        figure(ii);clf;
        
    
        surf(dat,'edgecolor','none')
        view([0 90])
        imagesc(dat);set(gca,'Ydir','normal')
        hold on
        plot(mt_learned,1:length(plot_inds),'k-o','linewidth',4);hold off
        if ii~=5
            set(gca,'ColorScale','log')
        end
        ylim([1 length(plot_inds)])
        xlim([1 mts_res])
        
        if ii==5
            clim([0.5 1])
            colorbar('ticklabelinterpreter','latex')
        else
            clim([0.01 1])
            set(gca,'ColorScale','log')
            colorbar('ticklabelinterpreter','latex','XTick', 10.^(-2.5:0.5:0))
        end
        set(gca,'fontsize',16,'ticklabelinterpreter','latex')
        xlabel('$\sigma_{\phi f}$','interpreter','latex')
        ylabel('$z(0)$ index','interpreter','latex') 
        drawnow
        saveas(gcf,['~/Desktop/expnum',num2str(expnum),'_ep',num2str(ep),'_mQ',num2str(mQ),'_stat',num2str(ii),'_mt.png'])
    end
end
end
end
