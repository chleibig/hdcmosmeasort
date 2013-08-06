function SaveResults( filename, ROIs )
%SAVERESULTS save results as event list
%   Detailed explanation goes here


fid = fopen(strcat(filename,'.events'),'wt');

%Header:
fprintf(fid, 'filename\tplace\ttime\tamplitude\tboss_column\tboss_row\tsensors\n');

unit_id = 0;
for l = 1:length(ROIs)
    for k = 1:length(ROIs(l).units)
        unit_id = unit_id + 1;
        for i = 1:length(ROIs(l).units(k).time)
            fprintf(fid,strcat(filename,'.h5','\t',num2str(unit_id),'\t',...
                num2str(ROIs(l).units(k).time(i)),'\t',...
                num2str(ROIs(l).units(k).amplitude(i)),'\t',...
                num2str(ROIs(l).units(k).boss_col),'\t',...
                num2str(ROIs(l).units(k).boss_row),'\t',...
                num2str(0),'\n'));
        end
    end
end

fclose(fid);

end
