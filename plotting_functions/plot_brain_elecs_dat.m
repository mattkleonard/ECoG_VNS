function [] = plot_brain_elecs_dat(dat, brain_dir,varargin)
%% Function plot electrodes on the surface of the brain with coloring
% that is proportional to the relative values of dat
% variable inputs included to make it possible to call this plot for
% variaous subjects and hemispheres.
%
% inputs:
% dat - num_electrodes long array. Will only plot data for non-zero
% elements of this array
%
% brain_dir: directory containing relevant anatomy files and electrode
% positions, and brain meshes.
%
% Variable inputs:
%
% 1 - subject name (default 'EC131') in general this will need to be
% changed for each subject
% 
% 2 - hemisphere (default 'lh') possible options ('rh', 'lh', 'both')
% entering 'both' will plot both the right and left hemisphere recons.
%
% 3 - brain alpha (value between 0 and 1) determines the transparency of
% the brain surface in the plot - transparent brains help plot depth
% electrodes
%
% 4 - sort electrodes alphabetically by anatomy
% if true, electrode coordinates will be sorted alphabetically by anatomy
% since the VNS data has been sorted alphabetically thus far, this is set
% to a default True, however if may be helpful in the future to change this
%
%% Load stuff:

subj = 'EC131';
hem = 'lh';
if length(varargin)>0
    if ~isempty(varargin{1})
    subj = varargin{1};
    end
end

%% The hemisphere - default set to 'lh'
% if user enters 'both', it will plot both hemispheres
if length(varargin)>1
     if ~isempty(varargin{2})
        hem = varargin{2};
     end
end
plot_full_brain = strcmpi(hem,'both');


%% Set brain alpha level
brain_alpha = 0.5;
if length(varargin)>2
    if ~isempty(varargin{3})
        brain_alpha = varargin{3};
    end
end

sort_elecs = true;
if length(varargin)>3
    if ~isempty(varargin{4})
        sort_elecs = varargin{4};
    end
end



load([brain_dir filesep subj filesep 'elecs/clinical_elecs_all.mat']); % load Electrode Map


if sort_elecs
    [~,order] = sort(anatomy(1:size(elecmatrix,1),4));
    elecmatrix = elecmatrix(order,:);
end
%% Get Colormap:
negative_rgb = [0 0 0.8];

positive_rgb = [0.8 0 0];

[cmap, data_axis] = make_approx_cbar_var(dat, negative_rgb, positive_rgb);

 sphere_size = 1.65;
 [sphx,sphy,sphz] = sphere(60);
 sphx = sphx*sphere_size; sphy = sphy*sphere_size; sphz = sphere_size*sphz;
% 

%% Normalize Data
dat_raw = dat;
mask = (dat ~= 0); % is True for significant channels 

dat = dat/(max(abs(dat(mask))));




% % roi = {'superiorfrontal','caudalmiddlefrontal','parstriangularis','parsopercularis','superiortemporal','precentral','postcentral','superiorparietal','supramarginal'};
% % clrs = linspace(-1,0,length(roi));

%roi_label_flag = 0;
plot_elecs_flag = 1;
%elec_label_flag = 0;
plot_bottom_side_views = 0;
elecSize = 80;
offset = 5;




if strcmpi(hem,'lh')
    offset = offset * -1;
end

%ctab_fid = fopen([rootdir '/MRI/aparc.annot.ctab']);
%ctab = textscan(ctab_fid,'%d%s%d%d%d%d');
%roi_names = ctab{2};
% roi = roi_names(2:end);
% clrs = linspace(-1,1,length(roi));

% % % % % % % rndColorOrder = randperm(length(clrs));
% % % % % % % clrs = clrs(rndColorOrder);
% % % % % % % clrs(find(clrs == 0)) = 0.25;
% % % % % % % 
% % % % % % % fid = fopen([rootdir '/MRI/' hem '.aparc.annot.dpv']);
% % % % % % % vert = textscan(fid,'%d%d%d%d%d');
% % % % % % % vert = vert{5};

%cmap = cbrewer('seq','Reds',101);



fig = figure;
% axes('Position', [0.85 0.20, 0.05,0.70])
% imagesc(repmat(flipud(data_axis'),1,30))
% colormap(cmap)
% ax = gca;
% ax.XTick = [];
% ax.YAxisLocation = 'right';
% ax.YTick = [1, length(data_axis)];
% ax.YTickLabel = {'Linguistic', 'Transitional'};
% axis equal tight
% 
% axes('Position', [0.15, 0.15, 0.75 0.75])
if ~plot_full_brain
    load([brain_dir filesep subj filesep 'Meshes' filesep subj '_' hem '_pial.mat']); % load Brain Mesh
    ctmr_gauss_plot_alpha(cortex,[0 0 0],0,hem,brain_alpha);
else
    load([brain_dir filesep subj filesep 'Meshes' filesep subj '_lh_pial.mat']);
    ctmr_gauss_plot_alpha(cortex,[0 0 0],0,'lh',brain_alpha);
    hold on;
    load([brain_dir filesep subj filesep 'Meshes' filesep subj '_rh_pial.mat']);
    ctmr_gauss_plot_alpha(cortex,[0 0 0],0,'rh',brain_alpha);
    hem = 'lh';
end
% % if plot_bottom_side_views
% %     subplot(1,2,1)
% %     view(270, 0);
% %     ctmr_gauss_plot(cortex,[0 0 0],0,hem);
% %     subplot(1,2,2)
% %     view(180, -60);
% % end
    

%view([-100,20])

kids = get(gca,'Children');



        
        
        

% % % % if roi_label_flag
% % % %     for i = 1:length(roi)
% % % %         roi_idx = find(strcmpi(roi_names,roi{i}));
% % % %         
% % % %         roi_verts{i} = find(vert == roi_idx);
% % % %         kids(2).FaceVertexCData(roi_verts{i},:) = repmat(clrs(i),length(roi_verts{i}),1);
% % % %     end
% % % % end
if plot_elecs_flag
    for i = 1:size(elecmatrix,1)
%         scatter3(elecmatrix(i,1)+offset,elecmatrix(i,2),elecmatrix(i,3),...
%             elecSize,'y','filled');
         if mask(i)
%           scatter3(elecmatrix(i,1)+offset,elecmatrix(i,2),elecmatrix(i,3),...
%              elecSize,cmap(round(dat(i)*100) + 1,:),'filled');
% 
                 [~,c_ind] = min(abs(data_axis - dat(i)));
                                    scatter3(elecmatrix(i,1)+offset,elecmatrix(i,2),elecmatrix(i,3),...
                       elecSize, cmap(c_ind,:),'filled'); %   colorscaling3(i) colorscaling2(i)],'filled');

                  % transp = abs(dat(i)/max(dat));
                    %transp = abs(dat_frame(i)/max(dat(:)));
                    %surf(sphx+elecmatrix(i,1)+offset,sphy+elecmatrix(i,2),sphz+elecmatrix(i,3), ...
                    %    'FaceColor', cmap(c_ind,:),'FaceAlpha', transp , 'LineStyle', 'none','FaceLighting','none', 'BackFaceLighting', 'unlit')
                    
         else
               scatter3(elecmatrix(i,1)+offset,elecmatrix(i,2),elecmatrix(i,3),...
               elecSize, [0.0 0.0 0.0]); %   colorscaling3(i) colorscaling2(i)],'filled');

          end
             
        hold on;
% % % % % %         if elec_label_flag
% % % % % %             text(elecmatrix(i,1)+offset,elecmatrix(i,2),elecmatrix(i,3),...
% % % % % %                 num2str(i),'Color','y');
% % % % % %         end
    end
end

% create_movie_time_bar(100)