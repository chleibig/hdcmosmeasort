function [ROIs, OL] = roisegmentation(data, metaData, paramsRoi, show)
%roisegmentation segments data into regions of interest
%
% Input
% =====
%
% data - (N_row x N_col x N_samples) array
% metaData - struct with fields describing data
%       *.sensor_rows - vector of sensor row IDs
%       *.sensor_cols - vector of sensor col IDs
%       *.frameStartTimes - N_samples vector of frame start times in ms
%       *.sr - sampling rate in kHz
%       *.filename_events - only needed for method 'cog'       
%
% paramsRoi - struct with following fields:
%       *.method = 'tce' or 'cog'
%       *.horizon - frames to the past and to the future of detected activity
%                   is taken for the temporal ROI
%       fields needed, if *.method = 'tce':
%       *.thr_factor - threshold in multiples of std
%       *.n_rows, *.n_cols, *.n_frames  environment in data coordinates
%       fields needed, if *.method = 'cog':
%       *.minNoEvents - the minimum # of events on a CoG sensor to be valid
%                       as a ROI seed
%       *.horizon - see under 'tce'
%
% show - flag to control graphical output
%
%
% Output
% ======
%
% ROIs - struct array with each struct element containing the relevant ROI
%        information in its fields:
%        *.X - (N x T) array containing a single
%        rectangular (in physical space) ROI
%        *.sensor_rows, *.sensor_cols - list of sensor coordinates
%                                for rectangular ROI
%        *.N_mask - (N x 1) boolean vector specifying which sensors
%        of rectangular X belong to actual roi [IN DATACOORDINATES!!]
%        *.T_mask - (T x 1) boolean vector for temporal subset
%        selection
%
%  OL -  (N_ROI x N_ROI) matrix with OL(k,l) being the fraction of
%        pixels being shared by the k-th and l-th region of interest
%
% christian.leibig@g-node.org, 16.07.13
%

switch paramsRoi.method
    case 'tce'
        %This method so far constructs a single ROI based on the presence
        %of (hyperspherical) threshold crossing events
        [ ROIs.X, ROIs.sensor_rows, ROIs.sensor_cols,...
            ROIs.T_mask, ROIs.N_mask ] = ...
            ROIIdentification(data, metaData.sensor_rows,...
            metaData.sensor_cols,paramsRoi, show);
        OL = [];
    case 'cog'
        [allRoi] = CoG_ROIs(metaData.filename_events,paramsRoi.minNoEvents,...
            max(metaData.sensor_rows), max(metaData.sensor_cols));
        %collect output
        for i = 1:allRoi.NumObjects
            %sensor row col pairs for each pixel in ROI_i:
            [sensor_rows_pxl_i,sensor_cols_pxl_i] = ...
                ind2sub(allRoi.ImageSize,allRoi.PixelIdxList{i});
            ROIs(i).sensor_rows = unique(sensor_rows_pxl_i);
            ROIs(i).sensor_cols = unique(sensor_cols_pxl_i);
            data_ROI_i = data(...
                (ROIs(i).sensor_rows(1) <= metaData.sensor_rows) & ...
                (metaData.sensor_rows <= ROIs(i).sensor_rows(end)),...
                (ROIs(i).sensor_cols(1) <= metaData.sensor_cols) & ...
                (metaData.sensor_cols <= ROIs(i).sensor_cols(end)),:);
            ROIs(i).X = reshape(data_ROI_i,...
                [size(data_ROI_i,1)*size(data_ROI_i,2) size(data_ROI_i,3)]);
            clear data_ROI_i
            %for *.N_mask we need the linear indices in data coordinates in
            %to ROIs(i).X
            [data_rows_pxl_i, data_cols_pxl_i] = sensorcoord2datacoord(...
                sensor_rows_pxl_i,sensor_cols_pxl_i,...
                ROIs(i).sensor_rows, ROIs(i).sensor_cols);
            lin_data_idx_pxl_i = sub2ind(...
                [length(ROIs(i).sensor_rows) length(ROIs(i).sensor_cols)],...
                data_rows_pxl_i, data_cols_pxl_i);
            ROIs(i).N_mask = ...
                false(length(ROIs(i).sensor_rows)*length(ROIs(i).sensor_cols),1);
            ROIs(i).N_mask(lin_data_idx_pxl_i) = true;
            %*.T_mask
            ROIs(i).time = allRoi.time{i};
            % it would be easy to convert time into frame indices, however,
            % it could be that time obtained from CoG_ROIs has a different
            % offset than the time used by the caller, hence check for that
            [ frameIdx ] = time2frame(ROIs(i).time,metaData.frameStartTimes);
            ROIs(i).T_mask = false(length(metaData.frameStartTimes),1);
            ROIs(i).T_mask(frameIdx) = true;
            ROIs(i).T_mask = fillvector( ROIs(i).T_mask, paramsRoi.horizon);
        end
        %overlap between rois:
        [ OL ] = roioverlap( allRoi.PixelIdxList );
        if show
            figure;colormap('gray');
            set(gcf,'position',get(0,'ScreenSize'));
            BackgroundAxes = axes('visible', 'off', 'units', 'normalized', 'Position', [0,0,1,1]);
            text(0.5, 0.99, 'Regions of interest', 'Parent', BackgroundAxes , ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','top');
            pltsize = ceil(sqrt(allRoi.NumObjects));
            set(gca,'Xtick',1:allRoi.NumObjects, 'Ytick',1:allRoi.NumObjects);
            for i = 1:allRoi.NumObjects
                subplot(pltsize,pltsize,i);
                imagesc(reshape(ROIs(i).N_mask,...
                    [length(ROIs(i).sensor_rows) length(ROIs(i).sensor_cols)]));
                set(gca,'XTickLabel',ROIs(i).sensor_cols(get(gca,'XTick')));
                xlabel('sensor cols');
                set(gca,'YTickLabel',ROIs(i).sensor_rows(get(gca,'YTick')));
                ylabel('sensor rows');
                title(strcat(num2str(i),': ',num2str(length(ROIs(i).time)),' events'));
            end
            
            figure;
            imagesc(OL);
            title(strcat('fraction of shared sensors (',...
                num2str(nnz(triu(OL))),' pairs > 0)'));
            xlabel('ROI index');
            ylabel('ROI index');
            axis square;
            colorbar;
        end
        
        
        
    otherwise
        error('Unknown method for region of interest segmentation.\n');
end



end