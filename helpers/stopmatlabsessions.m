function [ status, pids ] = stopmatlabsessions( multicorDir )
%[ status ] = stopmatlabsessions( multicorDir ) kills all slave processes
%by their pids save at startup time in multicorDir/pids.txt

% author: christian.leibig@g-node.org

pids = load([multicorDir filesep 'pids.txt'])';
fprintf('Shutting down %g matlab sessions.\n',length(pids));

delete([multicorDir filesep 'pids.txt']);
delete('setupslave.m');



if isunix
    [status, result] = system(['kill ' num2str(pids)]);    
end

if ispc
    [status, result] = system(['Taskkill /PID ' num2str(pids) ' /F']);
end

if ismac
    error('Mac compatibility not implemented yet.');
end

end