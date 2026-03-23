function [timet,errt] = gentab(datapath,writepath,pname,ptype)
%GENTAB Generates experiment tables
%   [TIMET,ERRT] = GENTAB(DATAPATH,WRITEPATH,PNAME,PTYPE)
%   Returns the experiment times and errors in TIMET and ERRT. 
%   PNAME is the problem e.g., 'vortex' and PTYPE is the type, e.g., 
%   'aligned'. Stores the tables such as miranda_times_alinged.txt 
%   with data from DATAPATH in WRITEPATH
%--------------------------------------------------------------------------
% 11/04/25, J.B., initial version
% 11/06/25, J.B., preparation for release
% 01/12/26, J.B., updated version
% 03/23/26, J.B., skip generation of 

% Add data to the path
addpath(datapath);

% probnames
probname = pname;

nstem_  = [datapath,'/',probname];
smps    = 50000;
if strcmp(ptype,'aligned')
    snames = {'pcg','direct_decoupled','direct'};
else
    snames = {'pcg','direct_nonsym'}; % ,'direct_sym'
end

nsolv   = length(snames);
tables  = cell(nsolv,1);

for i=1:nsolv
    ndata = [nstem_,'_',ptype,'_',snames{i},'_best.txt'];
    if ~exist(ndata,'file')      
        error(['gentable:: Data not found:',ndata,'. Please rerun',[' run_',ptype,'_',pname,'.m']]);
    end
    tables{i} = readtable(ndata,'NumHeaderLines',2);
end

% Data: "num runs" with 8 columns 
tbl     = cell(length(smps)*size(tables{1},1),1);
rnks    = tables{1}.rank;
nrnks   = length(rnks);

% tbl escape
tc = '&';
te = '\\';

timet       = zeros(nrnks,nsolv+1);
timet(:,1)  = rnks;
errt        = timet;
rowe        = zeros(1,nsolv);
rowt        = zeros(1,nsolv);

% loop through data
for i=1:length(smps)

    p = smps(i);


    cnt = 1;
    for ii = 1:nrnks % 1

        for jj=1:nsolv

            rowe(jj) = tables{jj}.rerr(ii);
            rowt(jj)= tables{jj}.time(ii);

        end

        timet(ii,2:end) = rowt;
        errt(ii,2:end)  = rowe;

        % format row
        [~,mix] = min(rowt);
        fl      = cell(nsolv,1);
        for iii=1:nsolv; fl{iii} = '{{'; end

        % annotate fastest row
        fc      = '}}';
        sbl     = '{\textbf{ \color{blue} ';
        fl{mix} = sbl;

        % multirow
        t1 = '';
        if ii==(1) && (strcmp(ptype,'unaligned') )
            t1 = ['\multirow{',num2str(nrnks),'}{*}{',num2str(p),'}',tc];
        end

        r = rnks(ii);
        
        tblr = [t1,num2str(r),tc];
        for iii=1:nsolv
            tblr = [tblr,fl{iii},num2str(rowt(iii),'%0.2f'),fc,tc]; %#ok<*AGROW>
        end
        for iii=1:nsolv
            if iii==nsolv; se=te; else; se=tc; end
            tblr = [tblr,num2str(rowe(iii),'%0.3f'),se];
        end
        %tblr = [tblr,te];
        
        % add hline
        if ii==(1)
            tblr = ['\hline ',tblr]; %#ok<AGROW>
        end

        tbl{(i-1)*(nrnks)+cnt} = tblr;

        cnt = cnt + 1;

    end

end

% save as text file
%writecell(tbl,[writepath,'/',[pname,'_table','_',ptype]]);

fname = fullfile([writepath,'/',[pname,'_times','_',ptype],'.txt']);
fid=fopen(fname,'w+');
fprintf(fid,'%s \n',sprintf('%% file: %s',mfilename));
fprintf(fid,'%s \n',sprintf('%% date: %s',string(datetime('now'))));
writetable(array2table(timet,'VariableNames',['rank',snames]),[writepath,'/',[pname,'_times','_',ptype],'.txt'],...
    'WriteMode','append','WriteVariableNames',true);
fclose(fid);

fname = fullfile([writepath,'/',[pname,'_errs','_',ptype],'.txt']);
fid=fopen(fname,'w+');
fprintf(fid,'%s \n',sprintf('%% file: %s',mfilename));
fprintf(fid,'%s \n',sprintf('%% date: %s',string(datetime('now'))));
writetable(array2table(errt,'VariableNames',['rank',snames]),[writepath,'/',[pname,'_errs' ,'_',ptype],'.txt'],...
    'WriteMode','append','WriteVariableNames',true);
fclose(fid);

end