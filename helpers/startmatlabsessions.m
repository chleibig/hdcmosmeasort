function startmatlabsessions( N, multicorDir)
%STARTMATLABSESSIONS Summary of this function goes here
%   Detailed explanation goes here

%generate m-file for parameter settings of independent matlab sessions to
%start (in case of directly passing variables into the system call, the 
%matlab sessions reside inside the starting session instead of being 
%independent)

% author: christian.leibig@g-node.org

if N <= 0; return; end

if isunix
    fprintf('Opening %g independent matlab sessions.\n',N);
    fileID = fopen('setupslave.m','w');
    fprintf(fileID,'function [pid] = setupslave()\n');
    fprintf(fileID,['addpath ' path ';\n']);
    fprintf(fileID,['fileID = fopen(' char(39) multicorDir filesep 'pids.txt' char(39) ',''a'');\n']);
    fprintf(fileID,'pid = feature(''getpid'');\n');
    fprintf(fileID,'fprintf(fileID,strcat(num2str(pid),''\\n''));\n');
    fprintf(fileID,'fclose(fileID);\n');
    fprintf(fileID,['startmulticoreslave(' char(39) multicorDir char(39) ');\n']);
    fprintf(fileID,'end\n');
    fclose(fileID);
    
    i = 0;
    while i < N
        i = i + 1;
        system([matlabroot, filesep, 'bin', filesep, './matlab',...
            ' -nodisplay -nosplash -nojvm -r "setupslave" &']);
    end
end

if ispc
    fprintf('Opening %g independent matlab sessions.\n',N);
    fileID = fopen('setupslave.m','w');
    fprintf(fileID,'function [pid] = setupslave()\n');
    currentPath = strrep(path,'\','\\');
    currentPath = [char(39) currentPath char(39)];
    currentPath = strrep(currentPath,';',[char(39) ',' char(39)]);
    fprintf(fileID,['addpath(' currentPath ');\n']);
    fprintf(fileID,['fileID = fopen(' char(39) multicorDir '\\' 'pids.txt' char(39) ',''a'');\n']);
    fprintf(fileID,'pid = feature(''getpid'');\n');
    fprintf(fileID,'fprintf(fileID,strcat(num2str(pid),''\\n''));\n');
    fprintf(fileID,'fclose(fileID);\n');
    fprintf(fileID,['startmulticoreslave(' char(39) multicorDir char(39) ');\n']);
    fprintf(fileID,'keyboard\n');
    fprintf(fileID,'end\n');
    fclose(fileID);
    
    i = 0;
    while i < N
        i = i + 1;
        system('matlab -nosplash -nojvm -r "setupslave" &');
    end
end

if ismac
    error('Mac compatibility not implemented yet!\n');
end

end