function SaveResults( filename, units )
%SAVERESULTS save results as basic event list
%   Detailed explanation goes here


fid = fopen(strcat(filename,'.basic_events'),'wt');

%Header:
fprintf(fid, 'filename\tplace\ttime\tamplitude\tboss_column\tboss_row\tsensors\n');

for k = 1:length(units)
    for i = 1:length(units(k).time)
        fprintf(fid,strcat(filename,'.h5','\t',num2str(k),'\t',...
            num2str(units(k).time(i)),'\t',...
            num2str(units(k).amplitude(i)),'\t',...
            num2str(units(k).boss_col),'\t',...
            num2str(units(k).boss_row),'\t',...
            num2str(0),'\n'));
    end
end

fclose(fid);

end
