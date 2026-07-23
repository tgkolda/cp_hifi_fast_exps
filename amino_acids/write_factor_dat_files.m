function write_factor_dat_files(M_in,prefix)

M = ktensor(M_in.lambda, M_in.U);
M = normalize(M,3);

%%
for k = 1:3
    fname = sprintf('%s-%d.dat',prefix,k);
    fid = fopen(fname,'w');
    % comments
    fprintf(fid,'# ');
    fprintf(fid,'Max value = %0.2e', max(max(M.U{k})));
    fprintf(fid,', ');
    fprintf(fid,'Min value = %0.2e', min(min(M.U{k})));
    fprintf(fid,'\n');
    % column headers
    if isa(M_in, 'ktensor_hifi')
        fprintf(fid,'x ');
    else
        fprintf(fid,'idx ');
    end
    for j = 1:size(M.U{k},2)
        fprintf(fid,'f%d ', j);
    end
    fprintf(fid,'\n');
    % data
    for i = 1:size(M.U{k},1)
        if isa(M_in, 'ktensor_hifi')
            fprintf(fid,'%g ', M_in.xvals{k}(i));
        else
            fprintf(fid,'%d ', i);
        end
        for j = 1:size(M.U{k},2)
            fprintf(fid,'%g ', M.U{k}(i,j));
        end
        fprintf(fid,'\n');
    end
    fclose(fid);
end

