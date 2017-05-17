function [frames] = plot_brain_elecs_sandbox(dat)%,dat2,dat3)

% dat = randn(1,256);
% dat = (dat - min(dat)) ./ (max(dat) - min(dat));
%
%dat = wieght_SVM_fs/(max(wieght_SVM_fs));
%dat = max_r2/max(max_r2);
%dat = Peak_LoxWeights_Sum_Norm;
%dat = Peak_HS_Weights_Sum_Norm;
mask = (dat ~= 0); % is True for significant channels 
normalize_dat = false;
if normalize_dat
    dat = dat./(max(dat));
end

scaling = max([max(dat), abs(min(min(dat)))]);
dat = dat/scaling;
%dat(dat<0) = dat(dat<0)/abs(min(min(dat)));
%dat(dat>0) = dat(dat>0)/(max(max(dat)));

%  dat2 = dat2/max(dat2);
%  dat3 = dat3/max(dat3);
% 
%  colorscaling1 = log(((dat - min(dat))*(exp(1)-1)/max(dat)+1));
%  colorscaling2 = log(((dat2 - min(dat2))*(exp(1)-1)/max(dat2)+1));
% colorscaling1 = sqrt(dat);
%   colorscaling2 = sqrt(dat2);
%   colorscaling3 = sqrt(dat3);
%%

subj = 'CH';
hem = 'lh';

roi = {'superiorfrontal','caudalmiddlefrontal','parstriangularis','parsopercularis','superiortemporal','precentral','postcentral','superiorparietal','supramarginal'};
clrs = linspace(-1,0,length(roi));

roi_label_flag = 0;
plot_elecs_flag = 1;
elec_label_flag = 0;
elecSize = 80;
offset = 5;

rootdir = '/Users/changlab/Documents/changrepo/matlab/analysis/ASL';

%%

load([rootdir '/MRI/Meshes/' subj '_' hem '_pial.mat']);
load([rootdir '/MRI/elecs/hd_grid.mat']);

if strcmpi(hem,'lh')
    offset = offset * -1;
end

ctab_fid = fopen([rootdir '/MRI/aparc.annot.ctab']);
ctab = textscan(ctab_fid,'%d%s%d%d%d%d');
roi_names = ctab{2};
% roi = roi_names(2:end);
% clrs = linspace(-1,1,length(roi));

rndColorOrder = randperm(length(clrs));
clrs = clrs(rndColorOrder);
clrs(find(clrs == 0)) = 0.25;

fid = fopen([rootdir '/MRI/' hem '.aparc.annot.dpv']);
vert = textscan(fid,'%d%d%d%d%d');
vert = vert{5};

%cmap = cbrewer('seq','Reds',101);

%%
% figure;
% range_dat= [min(min(dat)), max(max(dat))];
% for k = 1:size(dat,2)
%     dat_frame = squeeze(dat(:,k));
%     dat_frame = reshape(dat_frame,16,16);
%     imagesc(dat_frame)
%     caxis(range_dat)
%     colorbar
%      title({['Difference Between Transitional and Lexical High Gamma at Time = ' num2str(-2000+10*(k-0.5)), ' (ms)'];'Positive Correspondes to a Dominant Linguistic HG, Negative to Transitional HG'})
%  frames_squres(k) = getframe;
% 
%         
% end
% movie(frames_squres)



fig = figure;
ctmr_gauss_plot(cortex,[0 0 0],0,hem);
%view([-100,20])

kids = get(gca,'Children');

frames(length(100:300)) = struct('cdata',[],'colormap',[]);

for k = 100:300
    dat_frame = squeeze(dat(:,k));
    is_trans = dat_frame<0;
    
    colorscaling1 = dat_frame;
    colorscaling1(colorscaling1<0) = 0;
    colorscaling1 = sqrt(colorscaling1);
    
    colorscaling2 = dat_frame;
    colorscaling2(colorscaling2>0) = 0;
        colorscaling2 = abs(colorscaling2);
        colorscaling2 = sqrt(colorscaling2);
        
        
        

if roi_label_flag
    for i = 1:length(roi)
        roi_idx = find(strcmpi(roi_names,roi{i}));
        
        roi_verts{i} = find(vert == roi_idx);
        kids(2).FaceVertexCData(roi_verts{i},:) = repmat(clrs(i),length(roi_verts{i}),1);
    end
end
if plot_elecs_flag
    for i = 1:size(elecmatrix,1)
%         scatter3(elecmatrix(i,1)+offset,elecmatrix(i,2),elecmatrix(i,3),...
%             elecSize,'y','filled');
         if mask(i)
%          scatter3(elecmatrix(i,1)+offset,elecmatrix(i,2),elecmatrix(i,3),...
%              elecSize,cmap(round(dat(i)*100) + 1,:),'filled');
% 
            if is_trans(i);
            scatter3(elecmatrix(i,1)+offset,elecmatrix(i,2),elecmatrix(i,3),...
              elecSize, [colorscaling2(i) 0 0],'filled');
            else
                            scatter3(elecmatrix(i,1)+offset,elecmatrix(i,2),elecmatrix(i,3),...
              elecSize, [0 0 colorscaling1(i)],'filled');
            end

%             scatter3(elecmatrix(i,1)+offset,elecmatrix(i,2),elecmatrix(i,3),...
%               elecSize, [colorscaling1(i) colorscaling3(i) colorscaling2(i)],'filled');

%          else
%             scatter3(elecmatrix(i,1)+offset,elecmatrix(i,2),elecmatrix(i,3),...
%              elecSize,cmap(round(dat(i)*100) + 1,:));
          end
             
        hold on;
        if elec_label_flag
            text(elecmatrix(i,1)+offset,elecmatrix(i,2),elecmatrix(i,3),...
                num2str(i),'Color','y');
        end
    end
end
% 
 title(['High Gamma of mean Linguistic ERP at Time = ' num2str(-2000+10*(k-0.5)), ' (ms)'])
 frames(k-99) = getframe(fig);
end
film = movie(frames);