function lib=get_trigprod_lib(nstates,trigs,L,polys,include_crosstrig,add_polytrig)
    tags = [];
    if ~isempty(trigs)
        % trigs_nz = trigs(trigs~=0);
        trigs_nz = trigs;%(trigs~=0);
        for n=1:nstates/2
            tag_temp = zeros(length(trigs),nstates);        
            tag_temp_nz = zeros(length(trigs_nz),nstates);
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
        % tags = tags(~all(tags==0,2),:);
        tags = unique(tags,'rows');
    end
    if add_polytrig == 1
        tags_3 = prodtags(get_tags(1,[],nstates),tags(sum(tags~=0,2)==1,:));
    else
        tags_3 = [];
    end
    custom_tags = get_tags(polys,[],nstates);
    tags = [mat2cell(tags,ones(size(tags,1),1),size(tags,2));tags_3;mat2cell(custom_tags,ones(size(custom_tags,1),1),size(custom_tags,2))];
    lib = library('tags',tags);
end