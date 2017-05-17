function [sub_plot_coords] = position_subplot_grid(num_plots, varargin)
%% This function tiles a grid of subplots in a configuration that
% attempts to make the most effective use of space while preserving a
% reasonable aspect ratio of each subplot
% 
% This is a preliminary attempt. This follows a basic algorithm where the
% tiling finds the nearest square number and then adds additional rows and columns
% onto an n x n square, adding an extra column first and then adding an
% extra row.
%
% Inputs:
% num_plots - the number of subplots required
%
% variable inputs:
% tile_spacing - padding factor for the margins fraction of column/row size
% added as a buffer
spacing = 0.1;
if length(varargin)>0;
    if ~isempty(varargin{1})
        spacing = varargin{1};
    end
end


nearest_square = floor(sqrt(num_plots));
if nearest_square^2 == num_plots
    n_rows = nearest_square;
    n_cols = nearest_square;
else
    n_cols = nearest_square+1;
    n_rows = nearest_square;
    if n_cols*n_rows < num_plots
        n_rows = nearest_square+1;
    end
end

% List coords as a raster over columns then rows (1,2), (1,3),...(1,n),
% (2,1), ...(2,n),... (m,1),...(m,n);


% Set column/row tile sizes/spacing
% for n tiles, with n+2 margins of spacing = spacing*lengthtiles
size_col = 1/(n_cols+spacing*(2+n_cols)); % size of tiles along row direction
size_row = 1/(n_rows+spacing*(2+n_rows)); % size of tiles along col direction

margin_col = size_col*spacing;
margin_row = size_row*spacing;

col_origin = margin_col;
col_end = 1-(margin_col+size_col);
col_coords = linspace(col_origin,col_end, n_cols); % list of bottom corner coords of colums

row_origin = 1-(margin_row+size_row);
row_end = margin_row;
row_coords = linspace(row_origin, row_end,n_rows); % list of bottom corner row_coordinga

sub_plot_coords = zeros(n_cols*n_rows,4);
plot_ind = 1;
%figure;
for i = 1:n_rows
    for j = 1:n_cols
        sub_plot_coords(plot_ind,:) = [col_coords(j), row_coords(i), size_col, size_row];
        
 %       subplot('position', sub_plot_coords(plot_ind,:))
        plot_ind = plot_ind + 1;
        
    end
end


end
