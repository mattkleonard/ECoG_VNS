
%% 1) run segment anatomical MR in SPM5
% startup SPM (spm functions are used)
% coregister + reslice CT to anatomical MR
% in SPM5 reference image = MR, source image = CT
% segment MR
%% 2) generate surface for electrodes

% get_mask(subject,gray,white,outputdir,degree of smoothing,threshold) 
% e.g. get_mask(6,0.1) or get_mask(16,0.3)

% get_mask_V2('subjectX',... % subject name
%     './data/t1_class_2gray.nii',... % gray matter file
%     './data/t1_aligned_class.nii',... % white matter file
%     './data/t1_aligned.nii',... % t1 file
%     './data/',... % where you want to safe the file
%     'l',... % 'l' for left 'r' for right
%     13,0.2); % settings for smoothing and threshold

get_mask_V2('subjectX',... % subject name
    './data/t1_class_2gray.nii',... % gray matter file
    './data/t1_aligned_class.nii',... % white matter file
    './data/t1_aligned.nii',... % t1 file
    './data/',... % where you want to safe the file
    'l',... % 'l' for left 'r' for right
    17,0.2); % settings for smoothing and threshold

%% 3) select electrodes from ct
ctmr
% view result
% save image

%% 4) sort unprojected electrodes
sortElectrodes2;
% saves as electrodes_locx;

%% 5) plot electrodes 2 surface
% electrodes2surf(subject,localnorm index,do not project electrodes closer than 3 mm to surface)

% electrodes2surf(
    % 1: subject
    % 2: number of electrodes local norm for projection (0 = whole grid)
    % 3: 0 = project all electrodes, 1 = only project electrodes > 3 mm
    %    from surface, 2 = only map to closest point (for strips)
    % 4: electrode numbers
    % 5: (optional) electrode matrix.mat (if not given, SPM_select popup)
    % 6: (optional) surface.img (if not given, SPM_select popup)
    % 7: (optional) mr.img for same image space with electrode
    %    positions
% saves automatically a matrix with projected electrode positions and an image
% with projected electrodes
% saves as electrodes_onsurface_filenumber_inputnr2

% 1: for a grid use:
[out_els,out_els_ind]=electrodes2surf('name',...
    5,1,... % use these settings for the grid
    1:62,... % electrode numbers from the following file
    './data/electrodes_loc1.mat',... % file with electrode XYZ coordinates
    './data/subjectX_surface1_13_02.img',... % surface to which the electrodes are projected
    './data/t1_aligned.nii');
% 2: for a 2xN strip use:
[out_els,out_els_ind]=electrodes2surf(subject,4,1,1:16,'./data/electrodes_loc1.mat','./data/subjectX_surface1_13_02.img','./data/t1_aligned.nii');
% 3: for a 1xN strip use:
[out_els,out_els_ind]=electrodes2surf(subject,0,2,1:8,'./data/electrodes_loc1.mat','./data/subjectX_surface1_13_02.img','./data/t1_aligned.nii');



%% 6) combine electrode files into one and make an image

subject='name';
elecmatrix=nan(128,3);
switch subject
    case 'name'
    load('./data/name_electrodesOnsurface1_5.mat'); % F
    elecmatrix(1:62,:)=out_els;
%     load('./data/name_electrodesOnsurface1_4.mat'); % P
%     elecmatrix(49:64,:)=out_els;
%     load('./data/name_electrodesOnsurface2_4.mat'); % P
%     elecmatrix(65:80,:)=out_els;
  
    save('./data/name_electrodes_surface_loc_all1.mat','elecmatrix');
    [output,els,els_ind,outputStruct]=position2reslicedImage2(elecmatrix,'./data/t1_aligned.nii');

    for filenummer=1:100
        outputStruct.fname=['./data/electrodes_surface_all' int2str(filenummer) '.img' ];
        if ~exist(outputStruct.fname,'file')>0
            disp(['saving ' outputStruct.fname]);
            % save the data
            spm_write_vol(outputStruct,output);
            break
        end
    end
end

%% 6) generate cortex to render images:
gen_cortex_click_V2('name',0.1,[7 3],'l'); 


%% plot on surface

% load cortex
% load electrodes on surface
subject='name';
% load(['./data/name_cortex_sm7_3.mat']);
load(['./data/name_cortex.mat']);
ctmr_gauss_plot(cortex,[0 0 0],0)

load(['./data/name_electrodes_surface_loc_all1.mat']);
els=elecmatrix;
v_d=[270,0]; % view direction makes sure electrodes are visible on surface
a_offset=.1*max(abs(els(:,1)))*[cosd(v_d(1)-90)*cosd(v_d(2)) sind(v_d(1)-90)*cosd(v_d(2)) sind(v_d(2))];
els=els+repmat(a_offset,size(els,1),1);

loc_view(v_d(1),v_d(2))
el_add(els,'g',30);
% label_add(els);

% print('-painters','-r300','-dpng',strcat(['./figures/name_rot270_30_labels']));

% plot electrodes with sizes:
% use el_add_sizable
%-------------------------------------------------------------------------%


%%Brett Edits For Viewing
%all
el_add(els,'k',11);

%mesial
el_add(els(62:99,:),'k',8);
    
%mesial functional
el_add(els(62:99,:),'k',8);
el_add(els([86 94 95 96],:),'b',8);
el_add(els([74 92 93 94 95 96],:),'g',8);

% el_add(els(33:42,:),'b',8); %All PMC elecs
% el_add(els(39,:),'g',8);    %Episodic elecs        
% %el_add(els(33,:),'r',8);    %Rest elecs
% el_add(els([35 36 37 41 42],:),'k',8); %Bad chans

label_add(els);

%loc_light(90,0) %change numbers to address lighting

print('-painters','-r300','-dpng',strcat(['./data_SRb/figs/SRb_mesial_episodic']));

% plot electrodes with sizes:
% use el_add_sizable

background = get(gcf, 'color'); 
    % specify transparent background
    set(gcf,'color','none'); 
    % create output file
    set(gcf,'InvertHardCopy','off'); 
    print('-painters','-r300','-dpng', 'notTransparent.png');
    % read image data back in
    cdata = imread('notTransparent.png');
    % write it back out - setting transparency info
imwrite(cdata, 'SRb_GC_network_pairs.png', 'png', 'BitDepth', 16, 'transparency', background)

%----
%%New CTMR_Stanford Variations
% plot projected electrodes:
ctmr_gauss_plot(cortex,[0 0 0],0)
% scale to absmax
el_add_sizable(elecmatrix,[1:length(elecmatrix)]-round(length(elecmatrix)/2));
% or scale to max indicated 
el_add_sizable(elecmatrix,[1:length(elecmatrix)]-round(length(elecmatrix)/2),40);

loc_view(90,0)



%% testtest
ctmr_med_plot(cortex,cortex.vert(:,1)<-5)
loc_view(90,0)
