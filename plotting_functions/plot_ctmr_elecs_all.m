function [elec_colors] = plot_ctmr_elecs_all(rootdir,subj,hemi, elec_num_flag,elec_facealpha, varargin)

% subj = 'EC131';
% hemi = 'lh';
% elecmat = 'clinical';
% elec_num_flag = 1;
% facealpha = 0.1;
%
% Variable input;
% 1 - electrode ordering
permute_order = false;
if length(varargin) > 0
    if ~isempty(varargin{1})
        ordering = varargin{1};
        permute_order = true;
    end
end
% rootdir = '/Users/mattleonard/Documents/Research/dura/data_store2/imaging/subjects';

%%

%load([rootdir '/' subj '/Meshes/' subj '_' hemi '_pial.mat']);
%load([rootdir '/' subj '/elecs/' elecmat '_elecs_all.mat']);
load([rootdir '/' subj '_' hemi '_pial.mat']);
load([rootdir '/' 'clinical_elecs_all.mat']);

if permute_order
    elecmatrix=elecmatrix(ordering,:);
    eleclabels = eleclabels(1:size(elecmatrix,1),:);    
    eleclabels = eleclabels(ordering,:);
    anatomy = anatomy(1:size(elecmatrix,1),:);
    anatomy = anatomy(ordering,:);
    anatomy = strrep(anatomy, '_', '-');
end


% Set brain alpha:
% if isempty(brain_alpha) | ~isnumeric(brain_alpha)
%     brain_alpha = 1;
% end
    
figure;
ctmr_gauss_plot(cortex,[0 0 0],0,hemi);
hold on;

%array_names = regexp(eleclabels(:,1),'[0-9]','split');
%array_names = [array_names{:,1}];
array_names = anatomy(:,4);
[array_names_unique,idx] = unique(array_names(~strcmpi(array_names,'')),'stable');
array_names = array_names(~strcmpi(array_names,''));
elec_nums = regexp(eleclabels(:,1),['\d'],'match');

if max(idx) > size(elecmatrix,1)
    idx(find(idx > size(elecmatrix,1))) = [];
end
clrs = cbrewer('qual','Paired',length(array_names_unique));
elec_colors = zeros(size(eleclabels,1),3);
for i = 1:size(elecmatrix,1)
    %elec_array = regexp(eleclabels(i,1),'[0-9]','split');
    %elec_array = elec_array{:,1}(1);
    elec_array = anatomy{i,4};
    p(i) = scatter3(elecmatrix(i,1),elecmatrix(i,2),elecmatrix(i,3),50,clrs(find(strcmpi(elec_array,array_names_unique)),:),'filled');
    elec_colors(i,:) = clrs(find(strcmpi(elec_array,array_names_unique)),:);
    %p(i) = axes('position', [0.7 0.7 0.1 0.1])
    %plot(rand(100,1),randn(100,1),'r.')
    if elec_num_flag
        text(elecmatrix(i,1),elecmatrix(i,2),elecmatrix(i,3),char(elec_nums{i})')
    end
    hold on;
end
set(gcf,'Units','normalized');

kids = get(gca,'Children');
set(kids(end),'FaceAlpha',elec_facealpha);

legend(p(idx),array_names(idx),'FontSize',12,'Box','off');
% leg = axes('Position',[0.85 0.6 0.1 0.4]);
