%% Preproc final
% Make sure the functional files do not start with m
addpath('/mnt/d/toolboxes/spm12/spm12')

realign_done = 1;
tic
if realign_done == 0
clear;

addpath('/mnt/d/toolboxes/spm12/spm12/')
path = uigetdir('/mnt/d/codes/preprocessing related/');
cd (path);
Subj_list=dir('*0*');
nt = 404;
other_name_func = 1;
other_name_anat = 1;
com_dir = 0;
exclude_sub = zeros(1,length(Subj_list));
erorr_subjects = zeros(1,length(Subj_list));
errorText = cell(1,length(Subj_list));

tic
%%  Realignment
sub_not_processed = zeros(1,length(Subj_list));
method=questdlg('Do the want to process a subset of subjects comparing them to another directory?', 'Subset of Subjects','Yes','No','Yes');
switch method
        case 'Yes'
            path2 = uigetdir('/mnt/d/codes/preprocessing related/');
            other_dir = dir(path2);
            other_dir = other_dir([other_dir.isdir]);
            other_dir([1,2]) = [];   
            com_dir = 1;
        case 'No'
            exclude_sub = ones(1,length(Subj_list));   
end
% Exclude subjects not present in other directory
if com_dir == 1
    for subi =1:length(Subj_list)  
        for subj = 1:length(other_dir)
            if strcmp(Subj_list(subi).name,other_dir(subj).name)
                exclude_sub(subi) = 1;
                break                
            end
        end
    end
end

for subji =1:length(Subj_list)                                              % number of subjects

    if exclude_sub(subji) == 1                                              % Exclude subjects not present in other directory
        cd(fullfile(Subj_list(subji).folder,Subj_list(subji).name))             % Change directory to the subject folder
        Subj_dir = dir();                                                       % Subject folder
        Subj_dir([1,2]) = [];                                                   % Delete first two as they are not files
        if exist('func','dir')                                                  % if the func dir exits then make that as the functional folder
            func_fl_name = 'func';
            path_func = fullfile(Subj_dir(1).folder,'func');
        end
    
        if exist('anat','dir')                                                  % if the anat dir exits then make that as the anatomical folder
            anat_fl_name = 'anat';
            path_anat = fullfile(Subj_dir(1).folder,'anat');
        end
    
        if other_name_func && ~exist('func_fl_name','var')                         % if the path_func variable does not exits that means the name of functional folder is different
            % Specify the functional folder
            path_func = string(uigetdir(fullfile(Subj_list(subji).folder,Subj_list(subji).name),'Specify the Functional scan folder')); % Convert Char to str so we can use erase
            temp = erase(path_func,string(fullfile(Subj_list(subji).folder,Subj_list(subji).name)));             %Get the path of the func folder from the subject directory  %strsplit(path_func,'/');
            method=questdlg('Do the other subjects also have the same name', 'Same Name?','Yes','No','Yes');
            switch method
                case 'Yes'
                    func_fl_name = string(temp);
                    other_name_func = 0;
                case 'No'
                    other_name_func = 1;
                    clear func_fl_name
            end
        end    
        if other_name_anat && ~exist('anat_fl_name','var')                         % if the path_anat variable does not exits that means the name of anatomical folder is different
            % Specify the anatomical folder
            path_anat = string(uigetdir(fullfile(Subj_list(subji).folder,Subj_list(subji).name),'Specify the Anatomical scan folder'));
            temp = erase(path_anat,string(fullfile(Subj_list(subji).folder,Subj_list(subji).name))); %strsplit(path_anat,'/');
            method=questdlg('Do the other subjects also have the same name', 'Same Name?','Yes','No','Yes');
            switch method
                case 'Yes'
                    anat_fl_name = string(temp);
                    other_name_anat = 0;
                case 'No'
                    other_name_anat = 1;
                    clear anat_fl_name
            end
        end
        
        func_dir = dir(fullfile(Subj_list(subji).folder,Subj_list(subji).name,func_fl_name));                                              % Get the functional image folder
        func_dir([1,2]) = [];
        if exist('rest_name','var') 
            if strcmp(rest_name(end-2:end),'.gz')                          % if .gz file unzip it first
                rest_image = dir(fullfile(func_dir(1).folder,rest_name));
                gunzip(fullfile(func_dir(1).folder,rest_image.name)) 
                rest_image = dir(fullfile(func_dir(1).folder,rest_name(1:end-3)));
                temp = niftiread(fullfile(func_dir(1).folder,rest_image.name));
                nt = size(temp,4);
            elseif strcmp(rest_name(end-3:end),'.nii')
                rest_image = dir(fullfile(func_dir(1).folder,rest_name));
                temp = niftiread(fullfile(func_dir(1).folder,rest_image.name));
                nt = size(temp,4);
            end
            if isempty(rest_image)                                       % if the mentioned file does not exist
                sub_not_processed(subji) = 1;
                continue
            end
        else
            rest_image = dir(fullfile(func_dir(1).folder,'*.nii'));            
            if isempty(rest_image)
                rest_image = dir(fullfile(func_dir(1).folder,'*.gz'));
                if (length(rest_image) > 1)                                                                                             % If there are more than 1 gz files then
                    rest_image_name = uigetfile({'*.gz'},'More than one gz files in the Functional folder Please choose 1',path_func);  % ask for which gz file to choose
                    gunzip(fullfile(func_dir(1).folder,rest_image_name))                                                                   % Unzip the gz file
                    rest_image = dir(fullfile(func_dir(1).folder,rest_image_name(1:end-3)));                                               % Get the .nii file
                    temp = niftiread(fullfile(func_dir(1).folder,rest_image.name));
                    nt = size(temp,4);
                    method=questdlg('Do the other subjects also have the similar name', 'Similar Name?','Yes','No','Yes');
                    switch method
                    case 'Yes'
                        rest_name = char(inputdlg(['Please type the letters similar in other subjects. You selected ',rest_image_name]));
%                         other_name_rest = 0;
                    case 'No'
%                         other_name_rest = 1;
%                         clear rest_name
                    end
                end
            elseif (length(rest_image) > 1)                                                                                           % If there are more than 1 nifti files then
            clear rest_image
            rest_image_name = uigetfile({'*.nii'},'More than one nifti files in the Functional folder Please choose 1',path_func);
            rest_image = dir(fullfile(func_dir(1).folder,rest_image_name));                                               % Get the .nii file
                    method=questdlg('Do the other subjects also have the similar name', 'Similar Name?','Yes','No','Yes');
                    switch method
                    case 'Yes'
                        rest_name = inputdlg(['Please type the letters similar in other subjects. You selected ',rest_image_name]);
%                         other_name_rest = 0;
                    case 'No'
%                         other_name_rest = 1;
%                         clear rest_name
                    end
            end
            
        end
    
        for ii = 1:nt                                                           % number of timepoints
            pathlist{ii,1} = [char(fullfile(Subj_list(subji).folder,Subj_list(subji).name,func_fl_name,rest_image.name)),',',num2str(ii)];
        end
    
        matlabbatch{1}.spm.spatial.realign.estwrite.data = {pathlist}';
    
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;           % Quality 1--> max quality
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;                 % separation (in mm) between the points sampled in the reference image. Smaller speration gives more accurate results, but will be slower. 
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;                % FWHM Gaussian smothing kernel parameter in mm
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;                 % Register to mean
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;              % The method by which the images are sampled when estimating the optimum transformation. 2nd degree bspline
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];          % Directions in the volumes the values should wrap around in.
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';             % Optional weighting image to weight each voxel of the reference image differently when estimating the realignment parameters.
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];           % 
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;              % 4th degree bspline
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];          % no wraping
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;                % Because  of  subject  motion,  different images are likely to have different patterns of zeros from where it was not possible to sample data. With masking enabled, the program searches through the whole time series looking for voxels which need to be sampled from outside the original images. Where this occurs, that voxel is set to zero for the whole set of images (unless the image format can represent NaN, in which case NaNs are used where possible).
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';            % prefix
    
        spm('defaults', 'FMRI');
        spm_jobman('run', matlabbatch);
    %     cfg_util('run',matlabbatch);
        clear matlabbatch;
        clearvars -except path Subj_list nt other_name_func other_name_anat func_fl_name anat_fl_name path_func path_anat sc_st rest_name sub_not_processed com_dir exclude_sub erorr_subjects
    end 
end
	
%% Framewise Displacement 
% check the framewize displacement and exclude any if the criteria is not
% matched

excluded = 1 - exclude_sub;
p = 1;
fd_max_trans = zeros(length(Subj_list),1);
fd_max_rotat = zeros(length(Subj_list),1);

fd_mean_trans = zeros(length(Subj_list),1);
fd_mean_rotat = zeros(length(Subj_list),1);
for subji = 1:length(Subj_list)                                % For loop on number od subjects
    if exclude_sub(subji) == 1                                              % Exclude subjects not present in other directory
    
        % Get rest folder
        rest_folder = fullfile(Subj_list(subji).folder,Subj_list(subji).name,func_fl_name);            % folder path
        text_file = dir(fullfile(rest_folder,'*.txt'));
        rp_rest = load(fullfile(text_file.folder,text_file.name));              % input the text file generated at the Realignment stage
       
        rp_diff_trans = diff(rp_rest(:,1:3));                      % The first 3 parameters tell the displacement in x,y, and z direction in mm. Here the vector difference operator is used to get the derivative of vector. (framewise difference)
        rp_diff_rotat = diff(rp_rest(:,4:6)*180/pi);               % The last  3 parameters tell the rotation values pitch, yaw and roll in radians (here we convert them to degrees)
       
        multp = rp_diff_trans*rp_diff_trans';                      % diagonals will represent dx^2 + dy^2 + dz^2 
        fd_trans = sqrt(diag(multp));                              % sqrt (dx^2 + dy^2 + dz^2)
       
        multp = rp_diff_rotat*rp_diff_rotat';                      % diagonals will represent dpitch^2 + dyaw^2 + droll^2
        fd_rotat = sqrt(diag(multp));                              % sqrt (dpitch^2 + dyaw^2 + droll^2)
       
        fd_max_trans(subji,1) = max(fd_trans);                     % Get the maximum traslational framewize displacement
        fd_max_rotat(subji,1) = max(fd_rotat);                     % Get the maximum rotational framewize displacement
       
        fd_mean_trans(subji,1) = mean(fd_trans);                   % Get the mean traslational framewize displacement
        fd_mean_rotat(subji,1) = mean(fd_rotat);                   % Get the mean rotational framewize displacement
        
        if fd_max_trans(subji)>2 || fd_max_rotat(subji)>2 || fd_mean_trans(subji)>0.3 || fd_mean_rotat(subji)>0.3
            excluded(subji) = 1;
            excluded_subjects(p) = {Subj_list(subji).name};
            p = p+1;
        end
    end
end

% excluded = find(fd_max_trans>1.5|fd_max_rotat>1.5|fd_mean_trans>0.2|fd_mean_rotat>0.2);
% Try less stringent threshold for children data:
% excluded = find(fd_max_trans>2|fd_max_rotat>2|fd_mean_trans>0.3|fd_mean_rotat>0.3);
% excludedID = Subj_list(excluded);


if nnz(excluded)
    ask = questdlg(['A total of ',num2str(nnz(excluded)),' Subjects were discarded out of ',num2str(sum(exclude_sub)),' subjects. Do you want to Proceed'],'Proceed?','Yes','No','No');

    switch ask
        case 'No'
            return
    end   
end

end
%%
skull_strip=questdlg('Do you want to Skull Strip','Skull Strip', ...
                                                    'Yes AFNI','Yes FSL','No','No');

switch skull_strip 
    case 'Yes AFNI'
        sc_st = 2;         % 2: AFNI 
    case 'Yes FSL'
        sc_st = 1;         % 1: FSL
    case 'No'
        sc_st = 0;         % 0: No Skull Stripping
end

%% Preprocess the non-excluded subjects
subji_ex = 1;
T = zeros(length(Subj_list) - nnz(excluded),1);
for subji = 1:length(Subj_list)                                                       % For loop on number of subjects
    
    try
        if excluded(subji) == 0
            
            anat_folder = fullfile(Subj_list(subji).folder,...                                % Get the anatomical image folder
                            Subj_list(subji).name,anat_fl_name);
    %% Scull stripping
    % You may perform skull stripping at this stage for better coregisteration
    % output. For skull stripping you may either use AFNI or FSL
            if sc_st
                cd(anat_folder)                                                                   % Change Directory to the anat image folder
    
                if exist('anat_name','var')
                    if strcmp(anat_name(end-2:end),'.gz')                          % if .gz file unzip it first
                        anat_image = dir(fullfile(anat_name));
                        gunzip(full(anat_image.name))
                        anat_image = dir(fullfile(anat_name(1:end-3)));
                    elseif strcmp(anat_name(end-2:end),'.nii') 
                        anat_image = dir(fullfile(anat_name(1:end-3)));
                    end
                    if isempty(anat_image)                                       % if the mentioned file does not exist
                        sub_not_processed(subji) = 2;
                        continue
                    end
                    if sc_st == 1
                        system(['bet2 ',anat_image.name, [' b',anat_image.name]])                         % FSL command to skull strip the image
                    elseif sc_st == 2
                        system(['3dSkullStrip -input ',anat_image.name, ' -prefix' [' b',anat_image.name]]) % AFNI command to skull strip the image
                    end
                else
                    anat_image = dir('*.nii');                                                        % get the anatomical image
        
                    if isempty(anat_image)
                        anat_image = dir('*.gz');
                        
                        if (length(anat_image) > 1)                                                     % If there are more than 1 gz files
                            clear anat_image
                            anat_image_name = uigetfile({'*.gz'},...                                     % Ask the user to choose 1
                                'More than one nifti files in the anatomical folder Please choose 1',...
                                anat_folder);
                            gunzip(anat_image_name)                                                                   % Unzip the gz file
                            anat_image = dir(anat_image_name(1:end-3));                                               % Get the .nii file
                            method=questdlg('Do the other subjects also have the similar name', 'Similar Name?','Yes','No','Yes');
                            switch method
                            case 'Yes'
                                anat_name = char(inputdlg(['Please type the letters similar in other subjects. You selected ',anat_image_name]));
        %                         other_name_rest = 0;
                            case 'No'
        %                         other_name_rest = 1;
        %                         clear rest_name
                            end
                        elseif (length(anat_image) == 1) 
                            gunzip(anat_image.name)
                            anat_image = dir(anat_image.name(1:end-3));  
                        end
                        if sc_st == 1
                            system(['bet2 ',anat_image.name, [' b',anat_image.name]])                         % FSL command to skull strip the image
                        elseif sc_st == 2
                            system(['3dSkullStrip -input ',anat_image.name, ' -prefix' [' b',anat_image.name]]) % AFNI command to skull strip the image
                        end
        
                   elseif (length(anat_image) > 1)                                                     % If there are more than 1 nifti files
                            clear anat_image
                            anat_image_name = uigetfile({'*.nii'},...                                     % Ask the user to choose 1
                                'More than one nifti files in the anatomical folder Please choose 1',...
                                anat_folder);
                            anat_image = dir(anat_image_name);   
                        if sc_st == 1
                            system(['bet2 ',anat_image.name, [' b',anat_image.name]])                         % FSL command to skull strip the image
                        elseif sc_st == 2
                            system(['3dSkullStrip -input ',anat_image.name, ' -prefix' [' b',anat_image.name]]) % AFNI command to skull strip the image
                        end
                    end
                end
            end
    
    
    %% COREGISTRATION REST
    % Coregister anatomical image to the functional image
    
            func_folder = fullfile(Subj_list(subji).folder,...                                % Get the functional image folder
            Subj_list(subji).name,func_fl_name);
            cd(anat_folder)                                                                   % Change Directory to the anat image folder
            if sc_st == 1
                anat_image = dir('b*.nii.gz');     
                gunzip(fullfile(anat_image.folder,anat_image.name))
                anat_image = dir('b*.nii');     
                matlabbatch{1}.spm.spatial.coreg.estimate.ref = {[fullfile(anat_image.folder,anat_image.name),',1']};           % Reference image (this does not change) (Here its the anatomical image)
            elseif sc_st == 2
                anat_image = dir('b*.nii'); 
                matlabbatch{1}.spm.spatial.coreg.estimate.ref = {[fullfile(anat_image.folder,anat_image.name),',1']};       % Reference image (this does not change) (Here its the anatomical image)
    
            else
                anat_image = dir('*.nii'); 
                if length(anat_image) == 1
                    matlabbatch{1}.spm.spatial.coreg.estimate.ref = {[fullfile(anat_image.folder,anat_image.name),',1']};       % Reference image (this does not change) (Here its the anatomical image)
                else 
                    anat_image_name = uigetfile({'*.nii'},...                                     % Ask the user to choose 1
                        'More than one nifti files in the anatomical folder Please choose 1',...
                        anat_folder);
                    anat_image = dir(anat_image_name);
                    matlabbatch{1}.spm.spatial.coreg.estimate.ref = {[fullfile(anat_image.folder,anat_image.name),',1']};
                end
            end
            cd(func_folder)
            mean_func = dir('m*.nii');
            matlabbatch{1}.spm.spatial.coreg.estimate.source = {[fullfile(mean_func.folder,mean_func.name),',1']}; % Source image (This will change) (here its the functional image)
    
            func = dir('r*.nii');
            for ii= 1:nt    % Number of timepoints
                pathlist{ii,1} = strcat([fullfile(func.folder,func.name),',',num2str(ii)]);
            end
    
            matlabbatch{1}.spm.spatial.coreg.estimate.other = (pathlist);         % These are any images that need to remain in alignment with the moved image
    
            matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';  % Cost function is normalized mutual information
            matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];       % The average distance between sampled points (in mm).
            matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001]; % Iterations  stop  when differences between successive estimates are less than the required tolerance.
            matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [5 5];      % Histogram smoothing by Gaussian smoothing to apply to the 256x256 joint histogram. 
            
    %         cfg_util('run',matlabbatch);
            spm('defaults', 'FMRI');
            spm_jobman('run', matlabbatch);
            clear matlabbatch;
            
    %% Segmentation 
            matlabbatch{1}.spm.spatial.preproc.channel.vols = {fullfile(anat_image.folder,anat_image.name)};  % Select Volumes for processing
            matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
            matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
            matlabbatch{1}.spm.spatial.preproc.channel.write = [0 1];
            matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {'/mnt/d/toolboxes/spm12/spm12/tpm/TPM.nii,1'};
            matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
            matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {'/mnt/d/toolboxes/spm12/spm12/tpm/TPM.nii,2'};
            matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
            matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {'/mnt/d/toolboxes/spm12/spm12/tpm/TPM.nii,3'};
            matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
            matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {'/mnt/d/toolboxes/spm12/spm12/tpm/TPM.nii,4'};
            matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
            matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {'/mnt/d/toolboxes/spm12/spm12/tpm/TPM.nii,5'};
            matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
            matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {'/mnt/d/toolboxes/spm12/spm12/tpm/TPM.nii,6'};
            matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
            matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
            matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
            matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
            matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
            matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
            matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
            matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
            matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
            matlabbatch{1}.spm.spatial.preproc.warp.write = [0 1];
    %         cfg_util('run',matlabbatch);
            spm('defaults', 'FMRI');
            spm_jobman('run', matlabbatch);
            clear matlabbatch;
    
    %% Making CSF_MASK for REST
    
            func = dir('r*.nii');
            cd(anat_folder) 
            c3 = dir('c3*.nii');
            matlabbatch{1}.spm.util.imcalc.input = {
                [fullfile(func.folder,func.name),',1']%[path Subj_list(subji).name '/wrREST1.nii,1']
                 fullfile(c3.folder,c3.name) %[path Subj_list(subji).name '/wc3MPRAGE.nii']
                };
            matlabbatch{1}.spm.util.imcalc.output = 'CSF_MASK99.nii';
            matlabbatch{1}.spm.util.imcalc.outdir = {anat_folder};%{[path Subj_list(subji).name '/']};
            matlabbatch{1}.spm.util.imcalc.expression = 'i2>0.99';
            matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
            matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
            matlabbatch{1}.spm.util.imcalc.options.mask = 0;
            matlabbatch{1}.spm.util.imcalc.options.interp = 1;
            matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
    
    %         cfg_util('run',matlabbatch); 
            spm('defaults', 'FMRI');
            spm_jobman('run', matlabbatch);
            clear matlabbatch;
    
    
    %% Making WM_MASK for REST
    
            func_folder = fullfile(Subj_list(subji).folder,...                                % Get the functional image folder
                                  Subj_list(subji).name,func_fl_name);
            cd(func_folder)
            func = dir('r*.nii');
            temp = niftiread(fullfile(func.folder,func.name));
            nt = size(temp,4);
            anat_folder = fullfile(Subj_list(subji).folder,...                                % Get the anatomical image folder
                                  Subj_list(subji).name,anat_fl_name);
            cd(anat_folder) 
            c2 = dir('c2*.nii');
            matlabbatch{1}.spm.util.imcalc.input = {
            [fullfile(func.folder,func.name),',1']%[path Subj_list(subji).name '/wrREST1.nii,1']
             fullfile(c2.folder,c2.name) %[path Subj_list(subji).name '/wc3MPRAGE.nii']
            };
            matlabbatch{1}.spm.util.imcalc.output = 'WM_MASK99.nii';
            matlabbatch{1}.spm.util.imcalc.outdir = {anat_folder};%{[path Subj_list(subji).name '/']};
            matlabbatch{1}.spm.util.imcalc.expression = 'i2>0.99';
            matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
            matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
            matlabbatch{1}.spm.util.imcalc.options.mask = 0;
            matlabbatch{1}.spm.util.imcalc.options.interp = 1;
            matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
    
    %         cfg_util('run',matlabbatch);  
            spm('defaults', 'FMRI');
            spm_jobman('run', matlabbatch);
            clear matlabbatch;
    
    %% Extracting MEAN CSF and WM time-series
    
            cd(anat_folder) 
            % loading CSF MASKS
            CSF_MASK = spm_vol('CSF_MASK99.nii');%spm_vol([path Subj_list(subji).name '/CSF_MASK98.nii']);
            CSF_files  = spm_read_vols(CSF_MASK);
            image_dim = CSF_MASK.dim;
            CSF_Files_RS  = reshape(CSF_files,image_dim(1)*image_dim(2)*image_dim(3),1);
            clear image_dim
    
            % loading WM MASKS
            WM_MASK = spm_vol('WM_MASK99.nii');%spm_vol([path Subj_list(subji).name '/WM_MASK98.nii']);
            WM_files  = spm_read_vols(WM_MASK);
            image_dim = WM_MASK.dim;
            WM_Files_RS  = reshape(WM_files,image_dim(1)*image_dim(2)*image_dim(3),1);
    
            clear WM_REST CSF_REST REST_list REST_files image_dim2 REST_RS
            %%
    
        %     REST_folder = 'D:\NJIT\Research\Preprocessing\Practice\with_normalization\sub-10225\func';
            for imagei = 1:nt
    
                REST_list = spm_vol([fullfile(func.folder,func.name),',',num2str(imagei)]); %spm_vol([path Subj_list(subji).name '/wrREST1.nii,',num2str(imagei)]);
                REST_files = spm_read_vols(REST_list);
                image_dim2 = REST_list.dim;
                REST_RS = reshape(REST_files,image_dim2(1)*image_dim2(2)*image_dim2(3),1);
    
                WM_REST(:,imagei)  = REST_RS(find(WM_Files_RS > 0.5));
                CSF_REST(:,imagei) = REST_RS(find(CSF_Files_RS > 0.5));
                
    
            end
      %%    
            if isempty(WM_REST) || isempty(CSF_REST)
                erorr_subjects(subji) = 1;
                errorText{subji} = 'WM or CSF mask empty';
                continue
            end
            clear MEAN_WM_REST MEAN_CSF_REST VARname1 z*
    
            n1=isnan(CSF_REST);
            z1=find(n1>0);
            CSF_REST(z1)=0;
    
            n2=isnan(WM_REST);
            z2=find(n2>0);
            WM_REST(z2)=0;
    
            MEAN_WM  = mean(WM_REST);
            MEAN_CSF = mean(CSF_REST);
    
    
            MEAN_WM_REST  = MEAN_WM';
            MEAN_CSF_REST = MEAN_CSF';
    
            VARname1 = fullfile(func.folder,'covariance_wm_csf_REST.mat') ;%[path Subj_list(subji).name '/covariance_wm_csf_REST.mat'];
            save (VARname1, 'MEAN_WM_REST','MEAN_CSF_REST');
        
    
    %% TEMPORAL REGRESSION: removing nusisance signals and motion correction
    
            cd(func_folder)
    %         func = dir('r*.nii');
            % Loading One Image from MASKED_wrREST
            clear v_pair_* y_pair_* image_dim_* v_image_* y_image_* rp_* pca_CSF_* pca_WM_* mean_image_*
            v_pair_REST = spm_vol([fullfile(func.folder,func.name),',15']);%spm_vol([REST_folder 'wrREST1.nii,15']);
            y_pair_REST = spm_read_vols(v_pair_REST);
            image_dim_REST = v_pair_REST.dim;
            
		    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		    %%% *** Please Note: The REST_mask.nii is a mask to perform regression only within the brain ...
		    ... This script assumes that such a mask has already been created. One can use AFNI's 3dAutomask or ...
		    ... a modified version of FSL's bet module to do so on the functional data or any other suitable method that works ...
		    ... If the mask is not present in each subject's folder this part of the script will fail of course ...
		    ... Modify the script accordingly if you do not want to use any mask, which will obviously take much longer ....
		    
            system(['3dAutomask -prefix rest_mask.nii ',  func.name])
            % Loading REST_MASK
            MASK_list=spm_vol(fullfile(func.folder,'rest_mask.nii,1'));
            REST_MASK=spm_read_vols(MASK_list);
		    
		    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
            % Loading Time Series of MASKED_wrREST
            clear imagei
            for imagei = 1:nt
                v_image_REST = spm_vol(fullfile(func.folder,[func.name,',',num2str(imagei)]));%spm_vol([REST_folder 'wrREST1.nii,' num2str(imagei)]);
                y_image_REST(:,:,:,imagei) = spm_read_vols(v_image_REST);
            end
            
            fprintf('loading covariates for REST... \n')
            
            fprintf('First loading the motion parameters for REST... \n')
            txt_file = dir(fullfile(func.folder,'*.txt'));
            rp_REST1=load(fullfile(txt_file.folder,txt_file.name));
            rp_REST_temp=rp_REST1(1:nt,:);
            rp_REST = zscore(rp_REST_temp);
            rp_previous_REST = [0 0 0 0 0 0; rp_REST(1:end-1,:)];
            rp_auto_REST = [rp_REST rp_REST.^2 rp_previous_REST rp_previous_REST.^2];
            
            fprintf('Now loading csf and wm pca components for REST... \n')
            
            fprintf('Loading covariates from: \n')
            fprintf(fullfile(func.folder,'covariance_wm_csf_REST.mat'))
            cov_file=load(fullfile(func.folder,'covariance_wm_csf_REST.mat'));
            MEAN_WM_REST=cov_file.MEAN_WM_REST;
            MEAN_CSF_REST=cov_file.MEAN_CSF_REST;
            
            fprintf('Loading of covariates complete! \n')
            
            fprintf('regressing out covariates for REST... \n')
            mean_image_REST = mean(y_image_REST,4);   % calculate mean across time for all voxels | can be added back after regression to improve ICA performance
            
            %%
            clear vi vj vk
            for vi = 1:image_dim_REST(1)
                fprintf('%d',vi)
                for vj = 1:image_dim_REST(2)
                    for vk = 1:image_dim_REST(3)
                        if REST_MASK(vi,vj,vk)==0
                            
                            y_image_REST_regressed(vi,vj,vk,:) = zeros(1,1,1,nt);
                            
                        else
                            
                            % CSF, WM and MOTION Regressors and Regressors from OTHER SESSIONs for the same subject
                            
                            b = zscore([MEAN_WM_REST MEAN_CSF_REST rp_auto_REST]);
                            
                            [b,bint,y_image_REST_regressed(vi,vj,vk,:)] = regress(shiftdim(y_image_REST(vi,vj,vk,:),3),[b ones(nt,1)]);
                        end
                    end
                    
                end
            end
            
            fprintf('REST done! \n')
            fprintf('Release the Oompa Loompas!!\n')
            
            fprintf('saving REST image... \n')
            
            clear V
            V.dim   = v_pair_REST.dim;
            V.dt    = v_pair_REST.dt;
            V.mat   = v_pair_REST.mat;
            V.pinfo = v_pair_REST.pinfo;
            
            clear imagei
            for imagei = 1:nt
                V.n=[imagei 1];
                V.fname = fullfile(func.folder,['REG_',func.name]);
                V = spm_write_vol(V,y_image_REST_regressed(:,:,:,imagei)+mean_image_REST);
            end
            
    %% Nomralization
    
            matlabbatch{1}.spm.spatial.normalise.write.subj.def = {char(fullfile(anat_folder,['y_',anat_image.name]))};%{'/mnt/d/NJIT/Research/Preprocessing/Practice/new/sub-10471/anat/y_sub-10471_T1w.nii'};
            clear pathlist
            for ii= 1:nt    % Number of timepoints
                pathlist{ii,1} = strcat([fullfile(func.folder,['REG_',func.name]),',',num2str(ii)]);
            end
            matlabbatch{1}.spm.spatial.normalise.write.subj.resample = (pathlist);
            %%
            matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-90 -126 -72;  %[-78 -112 -70;78 76 85];
                                                                         90 90 108];
            matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
            matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
            matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
            
            cfg_util('run',matlabbatch); 
    %         spm('defaults', 'FMRI');
    %         spm_jobman('run', matlabbatch);
            clear matlabbatch
    %% Smoothing
    
            % Load Regressed Images
            for imagei = 1:nt
                reg_images{imagei, 1} = (fullfile(func.folder,['wREG_',func.name,',',num2str(imagei)]));
            end
    
            matlabbatch{1}.spm.spatial.smooth.data = reg_images;
            matlabbatch{1}.spm.spatial.smooth.fwhm = [6 6 6];
            matlabbatch{1}.spm.spatial.smooth.dtype = 0;
            matlabbatch{1}.spm.spatial.smooth.im = 0;
            matlabbatch{1}.spm.spatial.smooth.prefix = 's';
    
    %         cfg_util('run',matlabbatch);  
            spm('defaults', 'FMRI');
            spm_jobman('run', matlabbatch);
            clear matlabbatch
    
        end
        subji_ex = subji_ex + 1;
        disp('################################################################################################################')
        disp(['Completed Preprocessing Subject ',Subj_list(subji).name])
        T(subji_ex) = toc;
    
        disp(['Time Taken for subject ',Subj_list(subji).name, ' is ',num2str(T(subji_ex) - T(subji_ex - 1)),' seconds'])
        disp('################################################################################################################')
    catch exception
        erorr_subjects(subji) = 1;
        errorText{subji} = getReport(exception);
        continue
    end
end
    
    disp('Preprocessing is Finally over, I can now rest and you can get to work :-) ')
    
    %%
    % To format the participants.tsv file
    
    % input: file path---> the path of xlsx file
    
    
function T_formated = format_participants(filepath)

     T = readtable(filepath);
     T(:,6) = [];
     T_formated = unique(T,'rows');
end
