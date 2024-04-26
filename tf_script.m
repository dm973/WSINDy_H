% tf_rad = testfcn(Uobj(1),'meth','slDMD','param',0.5,'phifuns',phifun);
% tf1 = arrayfun(@(i)...
%         arrayfun(@(j)...
%             testfcn(Uobj(i),'stateind',j,'meth',tf_meth,'param',tf_param,'phifuns',dphi(phifun,sobol_inds(1)),...
%             'mtmax',floor(min(arrayfun(@(i)Uobj(i).dims-1,(1:length(Uobj))'))/8)),...
%         (1:nstates)','uni',0),(1:length(Uobj))','uni',0);
% 
arrayfun(@(U)U.getcorners,Uobj);
if and(isequal(tf_meth,'direct'),tf_param<0)
%     tf_param = max(arrayfun(@(j)max(arrayfun(@(i)get_tf_support(phifun,Uobj(j).dims,-tf_param,Uobj(j).ks(i)),1:nstates)),1:ntraj));
%     tf_param = get_tf_support(phifun,Uobj(1).dims,-tf_param,mean(arrayfun(@(j)1/mean(1./arrayfun(@(i)Uobj(j).ks(i),1:nstates)),1:ntraj)));
    tf_param = get_tf_support(phifun,Uobj(1).dims,-tf_param,mean(arrayfun(@(j)mean(arrayfun(@(i)Uobj(j).ks(i),1:nstates)),1:ntraj)));

end

tf1 = cell(size(Uobj));
for i=1:length(Uobj)
    ttemp_cell = cell(nstates,1);
    for j=1:nstates
        if mod(j,2)==1
            ttemp_cell{j} = testfcn(Uobj(i),'stateind',j,'meth',tf_meth,'param',tf_param,'phifuns',phifun,...
            'mtmax',floor(min(arrayfun(@(i)Uobj(i).dims-1,1:ntraj))*mtmax_fac/2),'subinds',subind);
        else
            ttemp_cell{j} = testfcn(Uobj(i),'stateind',j,'meth','direct','param',ttemp_cell{j-1}.rads,'phifuns',ttemp_cell{j-1}.phifuns,...
            'mtmax',floor(min(arrayfun(@(i)Uobj(i).dims-1,1:ntraj))*mtmax_fac/2),'subinds',subind);
        end
    end
    tf1{i} = ttemp_cell;
end

tf = [];
for k=1:length(sobol_inds)
    tf = [tf;arrayfun(@(i)...
            arrayfun(@(j)...
                testfcn(Uobj(i),'stateind',j,'meth','direct','param',tf1{i}{j}.rads,'phifuns',dphi(phifun,sobol_inds(k)),...
                'mtmax',floor(min(arrayfun(@(i)Uobj(i).dims-1,1:ntraj))/2),'subinds',subind),...
            (1:nstates)','uni',0),(1:ntraj)','uni',0)];
end
tf = [tf1;tf];
% arrayfun(@(i)arrayfun(@(j)disp(['rad (',num2str(i),',',num2str(j),')=',num2str(tf{i}{j}.rads)]),(1:nstates)','uni',0),(1:length(Uobj))','uni',0);
