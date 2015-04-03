function saveresults( ROIs, params )
% saveresults( ROIs, params ) 
% as event list, currently uses temporal offset from frameStartTimes!

filename = params.filenameResults;

tOffset = params.frameStartTimes(1);

fid = fopen(filename,'wt');

%Header:
fprintf(fid, ['filename\tplace\ttime\tamplitude\tboss_column\tboss_row'...
              '\tsensors\tquality\n']);

unit_id = 0;
for l = 1:length(ROIs)
    for k = 1:length(ROIs(l).units)
        unit_id = unit_id + 1;
        switch ROIs(l).units(k).state
            case 2
                quality = 'm';
            case 3
                quality = 's';
            otherwise
                quality = 'u';%for unchecked
        end
        for i = 1:length(ROIs(l).units(k).time) 
            fprintf(fid,strcat(filename,'\t',num2str(unit_id),'\t',...
                num2str(ROIs(l).units(k).time(i)+tOffset),'\t',...
                num2str(ROIs(l).units(k).amplitude(i)),'\t',...
                num2str(ROIs(l).units(k).boss_col),'\t',...
                num2str(ROIs(l).units(k).boss_row),'\t',...
                num2str(0),'\t',...
                quality,'\n'));
        end
    end
end

fclose(fid);

end
