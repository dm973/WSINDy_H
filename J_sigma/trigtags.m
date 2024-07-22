% freq = [[1,3];[3,1]];
% fmax = 4;
function tags = trigtags(freq,fmax)
    n = size(freq,1);
    nstates = size(freq,2);
    
    fs = -fmax:fmax;
    fs_c = repelem({fs},1,n);
    [fs_c{:}] = ndgrid(fs_c{:});
    fs_c = cell2mat(cellfun(@(f) f(:),fs_c,'un',0));
    tags = fs_c*freq;
    tags = unique(tags,'rows');
    
    C1 = arrayfun(@(i) tags(:,i)==-tags(:,i)',1:nstates,'un',0);
    C1 = all(cat(3,C1{:}),3);
    
    j=1;
    while j<length(C1)
        if any(C1(j,:))
            tags = tags(~C1(j,:),:);
            C1 = C1(~C1(j,:),~C1(j,:));
        end
        j = j+1;
    end
    tags = -tags;
    [~,I] = sort(vecnorm(tags,2,2),'ascend');
    tags = tags(I,:);
    tags = [(1+1i)*tags;(1-1i*tags(~all(tags==0,2),:))];
end