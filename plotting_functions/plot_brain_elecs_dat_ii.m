function [] = plot_brain_elecs_dat(dat, brain_dir)
%% Funcition plot electrodes on the surface of the brain with coloring
% that is proportional to the relative values of dat
%
%

%% Load stuff:

subj = 'EC131';
hem = 'lh';
load([brain_dir filesep subj '_' hem '_pial.mat']); % load Brain Mesh
load([brain_dir filesep 'clinical_elecs_all.mat']); % load Electrode Map

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
ctmr_gauss_plot(cortex,[0 0 0],0,hem);
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
               elecSize, [0 0 0]); %   colorscaling3(i) colorscaling2(i)],'filled');

          end
             
        hold on;
% % % % % %         if elec_label_flag
% % % % % %             text(elecmatrix(i,1)+offset,elecmatrix(i,2),elecmatrix(i,3),...
% % % % % %                 num2str(i),'Color','y');
% % % % % %         end
    end
end

% create_movie_time_bar(100)