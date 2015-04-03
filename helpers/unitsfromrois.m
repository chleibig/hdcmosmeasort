function [units] = unitsfromrois(ROIs)
%[units] = unitsfromrois(ROIs) adds ROI index for each unit in field 'k'
%
% christian.leibig@g-node.org, 2015-04-03
%

units = [];
for i = 1:length(ROIs);
    if ~isempty(ROIs(i).units);
        [ROIs(i).units.k] = deal(i);
        units = [units ROIs(i).units];
    end
end

end