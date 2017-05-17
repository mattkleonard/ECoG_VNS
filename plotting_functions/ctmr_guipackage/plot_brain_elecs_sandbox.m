function [] = plot_brain_elecs_sandbox(dat)



mask = (dat ~= 0); % is True for significant channels 

dat = dat/(max(abs(dat)));
colorscaling1 = dat;


subj = 'CH';
hem = 'lh';

roi = {'superiorfrontal','caudalmiddlefrontal','parstriangularis','parsopercularis','superiortemporal','precentral','postcentral','superiorparietal','supramarginal'};
clrs = linspace(-1,0,length(roi));

roi_label_flag = 0;
plot_elecs_flag = 1;
elec_label_flag = 0;
elecSize = 80;
offset = 5;


%%Changes!!!
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




fig = figure;
ctmr_gauss_plot(cortex,[0 0 0],0,hem);
%view([-100,20])

kids = get(gca,'Children');


        
        
        

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
%           scatter3(elecmatrix(i,1)+offset,elecmatrix(i,2),elecmatrix(i,3),...
%              elecSize,cmap(round(dat(i)*100) + 1,:),'filled');
% 
            if colorscaling1(i) > 0
            scatter3(elecmatrix(i,1)+offset,elecmatrix(i,2),elecmatrix(i,3),...
              elecSize, [colorscaling1(i) 0 0],'filled'); %   colorscaling3(i) colorscaling2(i)],'filled');
            else
                scatter3(elecmatrix(i,1)+offset,elecmatrix(i,2),elecmatrix(i,3),...
              elecSize, [0 0 -colorscaling1(i)],'filled'); %   colorscaling3(i) colorscaling2(i)],'filled');
            end
         else
              scatter3(elecmatrix(i,1)+offset,elecmatrix(i,2),elecmatrix(i,3),...
              elecSize, [0 0 0]); %   colorscaling3(i) colorscaling2(i)],'filled');

          end
             
        hold on;
        if elec_label_flag
            text(elecmatrix(i,1)+offset,elecmatrix(i,2),elecmatrix(i,3),...
                num2str(i),'Color','y');
        end
    end
end
% 
