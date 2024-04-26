disp(['-----getting library-----'])
disp(['libtype=',libtype,'.'])

if isequal(libtype,'polytrigprod')
%     L = 2;round(2*pi./mean(cell2mat(arrayfun(@(i)arrayfun(@(j)range(Uobj(i).Uobs{j}),1:nstates),(1:ntraj)','uni',0)),1));
    tags_1 = get_tags(polys,[],nstates);
    tags_2 = get_tags([],trigs'*min(L),nstates);
    tags_3 = prodtags(get_tags(1,[],nstates),tags_2);
    tags = [tags_1;tags_2];
    tags = mat2cell(tags,ones(size(tags,1),1),size(tags,2));
    tags = [tags;tags_3];
    lib = library('tags',tags);
elseif isequal(libtype,'trigprod')
%     L = 1;round(pi./mean(cell2mat(arrayfun(@(i)arrayfun(@(j)range(Uobj(i).Uobs{j}),1:nstates),(1:ntraj)','uni',0)),1));
    tags = [];
    if ~isempty(trigs)
        trigs_nz = trigs(trigs~=0);
        for n=1:Uobj(1).nstates/2
            tag_temp = zeros(length(trigs),Uobj(1).nstates);        
            tag_temp_nz = zeros(length(trigs_nz),Uobj(1).nstates);
            tags_1 = -1i*[tag_temp_nz(:,1:2*(n-1)) trigs_nz(:)*L(1) zeros(length(trigs_nz),1) tag_temp_nz(:,2*n+1:end)];
            tags_2 = -1i*[tag_temp_nz(:,1:2*(n-1)) zeros(length(trigs_nz),1) trigs_nz(:)*L(1) tag_temp_nz(:,2*n+1:end)];
            [ia, ib] = ndgrid(1:size(tags_1,1),1:size(tags_2,1));
            tags = [tags;tags_1(ia,:) + tags_2(ib,:)];
            tags_3 = 1i*[tag_temp(:,1:2*(n-1)) trigs(:)*L(1) zeros(length(trigs),1) tag_temp(:,2*n+1:end)];
            tags_4 = 1i*[tag_temp(:,1:2*(n-1)) zeros(length(trigs),1) trigs(:)*L(1) tag_temp(:,2*n+1:end)];
            [ia, ib] = ndgrid(1:size(tags_3,1),1:size(tags_4,1));
            tags = [tags;tags_3(ia,:) + tags_4(ib,:)];
            if include_crosstrig==1
                tags_5 = 1i*[tag_temp(:,1:2*(n-1)) trigs(:)*L(1) zeros(length(trigs),1) tag_temp(:,2*n+1:end)];
                tags_6 = -1i*[tag_temp_nz(:,1:2*(n-1)) zeros(length(trigs_nz),1) trigs_nz(:)*L(1) tag_temp_nz(:,2*n+1:end)];
                [ia, ib] = ndgrid(1:size(tags_5,1),1:size(tags_6,1));
                tags = [tags;tags_5(ia,:) + tags_6(ib,:)];
                tags_7 = -1i*[tag_temp_nz(:,1:2*(n-1)) trigs_nz(:)*L(1) zeros(length(trigs_nz),1) tag_temp_nz(:,2*n+1:end)];
                tags_8 = 1i*[tag_temp(:,1:2*(n-1)) zeros(length(trigs),1) trigs(:)*L(1) tag_temp(:,2*n+1:end)];
                [ia, ib] = ndgrid(1:size(tags_7,1),1:size(tags_8,1));
                tags = [tags;tags_7(ia,:) + tags_8(ib,:)];
            end
        end
        tags = unique(tags(~all(tags==0,2),:),'rows');
    end
    if add_polytrig == 1
        tags_3 = prodtags(get_tags(1,[],nstates),tags(sum(tags~=0,2)==1,:));
    else
        tags_3 = [];
    end
    custom_tags = get_tags(polys,[],nstates);
    tags = [mat2cell(tags,ones(size(tags,1),1),size(tags,2));tags_3;mat2cell(custom_tags,ones(size(custom_tags,1),1),size(custom_tags,2))];
    lib = library('tags',tags);
elseif isequal(libtype,'polytrigprod2')
    L = 1;round(pi./mean(cell2mat(arrayfun(@(i)arrayfun(@(j)range(Uobj(i).Uobs{j}),1:nstates),(1:ntraj)','uni',0)),1));
    fmax = max(trigs);
    [tags,~] = get_ndtrigtags(nstates,fmax,L,inf,0);
    tags = tags(sum(abs(tags),2)>0,:);
    tags_3 = prodtags(get_tags(1,[],nstates),tags);
    custom_tags = get_tags(polys,[],nstates);
    tags = [tags;custom_tags];
    tags = [mat2cell(tags,ones(size(tags,1),1),size(tags,2));tags_3];
    lib = library('tags',tags);
elseif isequal(libtype,'poly')
    tags = get_tags(polys,[],nstates);
    lib = library('tags',tags);
end

% lib.add_terms(term('fHandle',@(x1,x2,x3,x4) phi_I(sqr_fun(x1,x2,x3,x4))));
lib_length_io = length(lib.tags);
disp(['initial numterms=',num2str(lib_length_io)])

if toggle_trimlib == 1
    pt = mean(cell2mat(arrayfun(@(j)arrayfun(@(i)mean(Uobj(j).Uobs{i}),1:nstates),(1:ntraj)','uni',0)),1);
    beta = max(max(cell2mat(arrayfun(@(j)arrayfun(@(i)range(abs(Uobj(j).Uobs{i}-pt(i))),1:nstates),(1:ntraj)','uni',0))));
    tol = rangefun(beta);
    tic,
    lib.trimlib(pt,ord,tol,1,1);
    toc
end

tags = lib.tags;

Hlib = cellfun(@(tm) tm.fHandle, lib.terms, 'uni',0);
J = length(lib.terms);

disp(['numterms=',num2str(J)])