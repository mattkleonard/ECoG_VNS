
function [cmap, data_axis] = make_approx_cbar(dat)
%% Takes an array of values (dat) and creates a colormap that spans the range of the of values
% with negative values shaded blue, positive values shaded red with neutral
% values running to a light gray.
% the output is a cmap for a data axis. Points from dat should be plotted
% to the cmap index that contains the closest element of data_axis
dat_raw = dat;
pos_max_rgb = [56.5, 64.3, 21.6]/100;
neg_max_rgb = [44.3, 14.9, 43.1]/100;

%% Ignore values that are not relevant
mask_zeros = true;
if mask_zeros
    mask = (dat ~= 0); % is True for significant channels 
else
    mask = true(size(dat));
end

%% Normalize Data
maxval = (max(abs(dat(mask))))
dat = dat/maxval;

color_res = 500;    % Num of points in color gradation
gray_level = 0.5;   % Level of gray to fade into

data_axis = linspace(min(dat(mask)), max(dat(mask)), color_res);

if min(dat(mask)) >= 0  % set range from min to max all in red if all positive
    color_axis_red = linspace(gray_level, 1, color_res)';
    color_axis_green = linspace(gray_level,0, color_res)';
    color_axis_blue = linspace(gray_level,0,color_res)';
else %set range from min to max with <0 in blue
    red_pts = sum(data_axis >=0);
    blue_pts = sum(data_axis <0 );
    
    min_level = (1-gray_level)*abs(min(dat(mask)));
    max_level = (1-gray_level)*abs(max(dat(mask)));
    
    color_axis_red = [linspace(1 - gray_level - min_level, gray_level, blue_pts)'; linspace(gray_level, gray_level+max_level, red_pts)'];
    color_axis_green = [linspace(1 - gray_level - min_level, gray_level,blue_pts)'; linspace(gray_level, 1 - gray_level - max_level, red_pts)'];
    color_axis_blue = [linspace(gray_level+min_level,gray_level, blue_pts)'; linspace(gray_level, 1 - gray_level - max_level, red_pts)'];

    
%     color_axis_red = [linspace(0.8*abs(min(dat(mask))), gray_level, blue_pts)'; linspace(gray_level, 0.8*abs(max(dat(mask))), red_pts)'];
%     color_axis_green = [linspace(0.25, gray_level,blue_pts)'; linspace(gray_level, 0.8*abs(max(dat(mask))), red_pts)'];
%     color_axis_blue = [linspace(0.8*abs(min(dat(mask))), gray_level, blue_pts)'; linspace(gray_level, 0.25*abs(max(dat(mask))), red_pts)'];

%     color_axis_red = [linspace(neg_max_rgb(1), gray_level, blue_pts)'; linspace(gray_level, 0.565, red_pts)'];
%     color_axis_green = [linspace(neg_max_rgb(2), gray_level,blue_pts)'; linspace(gray_level, .64, red_pts)'];
%     color_axis_blue = [linspace(0.44*abs(min(dat(mask))), gray_level, blue_pts)'; linspace(gray_level, 0.21, red_pts)'];

end

cmap = [color_axis_red, color_axis_green, color_axis_blue];

%% Display the colormaping
display_cbar = true;
if display_cbar
    cvals = linspace(min(dat_raw(mask)), max(dat_raw(mask)), color_res); % used to display colorbar
    figure; imagesc(flipud(cvals'))
    colormap(cmap)
    colorbar
end

end

