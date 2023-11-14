function [all_TE, pfilePath, Trecon] = QSM_recon(bnumber, tnumber, output_path, pfilePath)

% =========================================================================
% Input
% --------------
% bnumber           : patient ID
% tnumber           : exam ID
%     	- default: parsed from pwd
% output_path       : output path for swan data
%       - default: bnumber/tnumber/swan_qsm
%
% Output
% --------------
% Trecon        	: reconstruction time
%
% Description: This is a wrapper function to run QSM reconstruction
%
% Jingwen Yao
% jingwen.yao@ucsf.ucsf.edu
% Date created: 8 July 2021
% Date modified: 8 July 2021
% Modified from Eason's code
% =========================================================================

%% add path

addpath('/netopt/share/lib/local/brain/matlab/');
addpath(genpath('/netopt/share/lib/local/brain/matlab/arc'));
addpath(genpath('/home/jyao3/030_QSM/01_Code'));

if nargin < 4
    
    %% set default
    
    path_root = '/data/7T_hunt';
    
    if nargin < 1
        currentPath = pwd;
        C = strsplit(currentPath,'/');
        if length(C) < 2; error('Wrong directory!'); end;
        tnumber = C{end};
        bnumber = C{end-1};
        if ~strcmp(bnumber(1),'b') && ~strcmp(bnumber(1:4),'temp'); error('Wrong directory!'); end;
        if ~strcmp(tnumber(1),'t') && ~strcmp(tnumber(1:3),'for'); error('Wrong directory!'); end;
        output_path = [path_root '/' bnumber '/' tnumber '/swan_qsm'];
    end
    
    if nargin < 3
        output_path = [path_root '/' bnumber '/' tnumber '/swan_qsm'];
    end
    
    fprintf('## Reconstructing patient %s exam %s \n', bnumber, tnumber);
    fprintf('## QSM will be saved in %s \n', output_path);
    
    mkdir(output_path);
    
    %% Recognize and read pfile
    
    pfilePath = [];
    examPath = [path_root '/' bnumber '/' tnumber];
    % go through everything within the folder or one layer beyond
    
    Efolder = dir(examPath);
    Efolder(~[Efolder.isdir]) = [];
    Efolder(~cellfun(@isempty,strfind({Efolder.name},'.'))) = [];
    % Efolder(cellfun(@isempty,strfind({Efolder.name},'E'))) = [];
    for jj = 1:length(Efolder)
        Epath = [examPath '/' Efolder(jj).name];
        seqFolder = dir(Epath);
        seqFolder(~[seqFolder.isdir]) = [];
        seqFolder(cellfun(@isempty,strfind({seqFolder.name},'raw'))) = [];
        for ii = 1:length(seqFolder)
            seqFile = dir([Epath '/' seqFolder(ii).name]);
            seqFile([seqFile.isdir]) = [];
            if isempty(seqFile)
                continue;
            end
            seqPath = [Epath '/' seqFolder(ii).name '/' seqFile(1).name];
            [~,outstr] = system(['printraw ' seqPath '| grep ''se_desc''']);
            if ~isempty(strfind(outstr, 'SWAN'))
                pfilePath = seqPath;
            end
        end
        if isempty(pfilePath)
            seqFile = dir(Epath);
            seqFile([seqFile.isdir]) = [];
            if isempty(seqFile)
                continue;
            end
        end
        if isempty(pfilePath)
            for ii = 1:length(seqFile)
                seqPath = [Epath '/' seqFile(ii).name];
                [~,outstr] = system(['printraw ' seqPath '| grep ''se_desc''']);
                if ~isempty(strfind(outstr, 'SWAN'))
                    pfilePath = seqPath;
                end
            end
        end
    end
    
    if isempty(pfilePath)
        seqFile = dir(examPath);
        seqFile([seqFile.isdir]) = [];
        for ii = 1:length(seqFile)
            seqPath = [examPath '/' seqFile(ii).name];
            [~,outstr] = system(['printraw ' seqPath '| grep ''se_desc''']);
            if ~isempty(strfind(outstr, 'SWAN'))
                pfilePath = seqPath;
            end
        end
    end
    
    if isempty(pfilePath)
        seqFile = dir([examPath '/../orig_raw']);
        seqFile([seqFile.isdir]) = [];
        for ii = 1:length(seqFile)
            seqPath = [examPath '/../orig_raw/' seqFile(ii).name];
            [~,outstr] = system(['printraw ' seqPath '| grep ''se_desc''']);
            if ~isempty(strfind(outstr, 'SWAN'))
                pfilePath = seqPath;
            end
        end
    end
    
    if isempty(pfilePath); error('Cannot find SWAN P file!'); end;
    
end

fprintf('## P file path: %s \n', pfilePath);

%% load pfile

tic;

hdr = read_raw_hdr_26x_Jingwen(pfilePath);
load_hdr_fields_fast;
yres = scan_phase;
voxel_size = [hdr.dfov/xres, hdr.dfov/yres, hdr.slthick];

nEcho = hdr.rdb_hdr_nechoes;
fprintf('--> Number of echoes: %i \n', nEcho);

all_TE = zeros(1,nEcho);
[~,outstr] = system(['printraw ' pfilePath '| grep ''rdb_hdr_echotimes''']);
C = strsplit(outstr);
for ii = 1:nEcho
    all_TE(ii) = str2double(C(ii*3))*10^3;
end
fprintf('--> Echo times: %.3f \n', all_TE);

[~,outstr] = system(['printraw ' pfilePath '| grep ''magstrength''']);
C = strsplit(outstr);
B0 = str2double(C(3))/10^4;
fprintf('--> B0: %.1f Tesla \n', B0);

%% Reconstruction

if 1 %isempty(dir(sprintf('%s/*_magni_echo%d.idf', output_path, nEcho))) ...
%         && isempty(dir(sprintf('%s/idf_echoes/*_magni_echo%d.idf', output_path, nEcho)))
    for echo = 1:nEcho
        % read k space data
        [~, k] = read_raw_data_swan(pfilePath, 1, echo);
        khat = k(:,:,:,first_slice:last_slice);
        % ARC recon
        khat = recon_arc_kdata_par(khat);
        % reconstruct image by ifft2c
        im = ifft2c(khat);
        % flip and exchange the dimensions for correct direction in idf file
        im = flip(im,4);
        im = permute(im,[2 1 4 3]); % row x col x slice x coil
        
        % combine coils
        if echo == 1
            magni1 = abs(im); phase1 = angle(im); clear im;
        else
            magni2 = abs(im); phase2 = angle(im); clear im;
            [mcpcs1, mcpcs2] = MCPC_3D_S(magni1, magni2, phase1, phase2, all_TE(1), all_TE(echo), voxel_size);
            unwrapped1 = laplacian_unwrap(mcpcs1, voxel_size, [0,0,6]);
            unwrapped2 = laplacian_unwrap(mcpcs2, voxel_size, [0,0,6]);
            
            % save magnitude and phase
            magni_output_filepath = fullfile(output_path, sprintf('%s_magni_echo%d', tnumber, 1));
            magni_output = create_idf_from_raw_hdr(hdr, magni_output_filepath, 3);
            magni_output.img = sqrt(sum(magni1.^2, 4));
            write_idf_image(magni_output_filepath, magni_output);
            
            phase_output_filepath = fullfile(output_path, sprintf('%s_phase_echo%d', tnumber, 1));
            phase_output = create_idf_from_raw_hdr(hdr, phase_output_filepath, 7);
            phase_output.img = mcpcs1;
            write_idf_image(phase_output_filepath, phase_output);
            
            unwrapped_output_filepath = fullfile(output_path, sprintf('%s_unwrap_echo%d', tnumber, 1));
            unwrapped_output = create_idf_from_raw_hdr(hdr, unwrapped_output_filepath, 7);
            unwrapped_output.img = unwrapped1;
            write_idf_image(unwrapped_output_filepath, unwrapped_output);
            
            % save magnitude and phase
            magni_output_filepath = fullfile(output_path, sprintf('%s_magni_echo%d', tnumber, echo));
            magni_output = create_idf_from_raw_hdr(hdr, magni_output_filepath, 3);
            magni_output.img = sqrt(sum(magni2.^2, 4));
            write_idf_image(magni_output_filepath, magni_output);
            
            phase_output_filepath = fullfile(output_path, sprintf('%s_phase_echo%d', tnumber, echo));
            phase_output = create_idf_from_raw_hdr(hdr, phase_output_filepath, 7);
            phase_output.img = mcpcs2;
            write_idf_image(phase_output_filepath, phase_output);
            
            unwrapped_output_filepath = fullfile(output_path, sprintf('%s_unwrap_echo%d', tnumber, echo));
            unwrapped_output = create_idf_from_raw_hdr(hdr, unwrapped_output_filepath, 7);
            unwrapped_output.img = unwrapped2;
            write_idf_image(unwrapped_output_filepath, unwrapped_output);
        end
    end
end

%% convert of nifti
if exist([output_path '/Phase.nii.gz'], 'file') ~=2
    
    cd(output_path);
    % put phase tissue in folder
    if ~isempty(dir(sprintf('%s/*_tissue_phase_echo%d.idf', output_path, nEcho)))
        mkdir([output_path '/idf_echoes/idf_tissue_phase']);
        cmd = sprintf('mv %s %s', ...
            [output_path '/*_tissue_phase_echo*.idf'], [output_path '/idf_echoes/idf_tissue_phase']);
        system(cmd);
    elseif ~isempty(dir(sprintf('%s/idf_echoes/*_tissue_phase_echo%d.idf', output_path, nEcho)))
        mkdir([output_path '/idf_echoes/idf_tissue_phase']);
        cmd = sprintf('mv %s %s', ...
            [output_path '/idf_echoes/*_tissue_phase_echo*.idf'], [output_path '/idf_echoes/idf_tissue_phase']);
        system(cmd);
    end
    % convert images to nifti
    if ~isempty(dir(sprintf('%s/*_magni_echo%d.idf', output_path, nEcho)))
        for nn = 1:nEcho
            cmd = sprintf('i2nii -o %s -z y %s', ...
                output_path, ['*_magni_echo' num2str(nn) '.idf']);
            system(cmd);
            cmd = sprintf('i2nii -o %s -z y %s', ...
                output_path, ['*_phase_echo' num2str(nn) '.idf']);
            system(cmd);
        end
    elseif ~isempty(dir(sprintf('%s/idf_echoes/*_magni_echo%d.idf', output_path, nEcho)))
        for nn = 1:nEcho
            cmd = sprintf('i2nii -o %s -z y idf_echoes/%s', ...
                output_path, ['*_magni_echo' num2str(nn) '.idf']);
            system(cmd);
            cmd = sprintf('i2nii -o %s -z y idf_echoes/%s', ...
                output_path, ['*_phase_echo' num2str(nn) '.idf']);
            system(cmd);
        end
    end
    
    % combine to one nifti
    system('fslmerge -t Magni.nii.gz *magni_echo*.nii.gz');
    system('fslmerge -t Phase.nii.gz *phase_echo*.nii.gz');
    mkdir('nii_echoes');
    system(sprintf('mv %s_*.nii.gz nii_echoes', tnumber));
end

%% end
Trecon = toc;
fprintf('Total time: %f \n', Trecon);

rmpath('/netopt/share/lib/local/brain/matlab/');

end

%% helper functions

function im = ifft2c(k)
im = ifftshift(ifft2(ifftshift(k)));
end

function [hdr, k] = read_raw_data_swan(pfile, fftflag, echo)

% Read_raw_data reads pfile from the GE scanner and assembles raw data to a
% K-Space matrix
%
% Input:  --pfile    name of raw file
%         --fftflag  flag to do IFFT along the slice direction (for 3D
%                    acquisition)
%         --necho    echo number if multiple echoes were acquired
%
% Wei Bian Oct, 2012 (Modified from the code given by Suchi)

% Xiaowei Zou, add ifftshift before applying fourier transform in slice direction to eliminate pi phase jump in slice direction

% Eason Chen, 08/30/2017, modified for GE Ax SWAN sequence with ARC

[fid,message] = fopen(pfile, 'r','l');

if ~isempty(message)
    display(['Raw data file can NOT be opened: ' message]);
    return
end

switch nargin
    case 1
        fftflag = 0; % Do Not apply IFFT along the slice direction if fftflag not assigned
        echo = 1;   % read raw data for the first echo if necho not assigned
    case 2
        echo = 1;
end

% read header
hdr=read_raw_hdr_26x_Jingwen(pfile);
load_hdr_fields_fast;
yres = scan_phase; % By Eason: the kspace data in pfile is already zero-filled.

baseline = xres*2*pt_size;
if pt_size == 2
    type_flag = 'int16';
else
    type_flag = 'int32';
end
k = zeros(xres,yres,ncoil,nlocslab);

% ncoil = 26;

% skip the header
fseek(fid,hdr_size,'bof');
% search the start byte for the required echo
fseek(fid, xres*2*(yres+1)*pt_size*(echo-1), 0); %
for nc=1:ncoil  % loop for coils
    fprintf('Loading coil %d', nc);
    for ns=1:nlocslab  % loop for slices
        %         fprintf(' %d', ns);
        % skip baseline view
        fseek(fid, baseline, 0);
        % read raw data
        raw = fread(fid, xres*2*yres, type_flag);
        complex_raw = raw(1:2:end) + 1i*raw(2:2:end);
        k(:,:,nc,ns) = reshape(complex_raw,xres,yres);
        fseek(fid, xres*2*(yres+1)*pt_size*(necho-1), 0);
    end
    fprintf('\n');
end

% do IFFT along slice direction
if fftflag
    for nc=1:ncoil
        k(:,:,nc,:) = ifft(ifftshift(squeeze(k(:,:,nc,:)),3),[],3);
    end
end

fclose(fid);

end


%%
function [combined1, combined2] = MCPC_3D_S(magni1, magni2, phase1, phase2, TE1, TE2, spatial_res)
% input:
%   cmplx1: the complex image of first echo. should be a 4D matrix:
%           [y X x X z X coil]
%   cmplx2: the comples image of second echo.
%   TE1, TE2: echo time of first and second echo
%   spatial_res: spatial resolution of the input images, in mm.
% output:
%   combined1, combined2: combined phase of echo1 and echo2
% Note: might need a lot of memory, especially running on 32-ch data

% addpath('../unwrap');
% addpath('../utils');

% calculate the phase difference between echoes
dP = angle(sum(magni2 .* magni1 .* exp(1i * (phase2-phase1)), 4)); % delta phase
% unwrap the phase difference
dPu = laplacian_unwrap(dP, spatial_res, [0,0,6]); % delta phase unwrapped
% high res estimate of coil phase offset
p0c = angle(exp(1i * (phase1 - TE1/(TE2-TE1) * repmat(dPu, 1,1,1,size(magni1, 4)))));
% smoothing the estimate of coil offset
p0cs = lp_filter(magni1, p0c, spatial_res); % smoothed estimate of coil offset

% combine phases using the smooth estimate of coil offset
combined1 = angle(sum(magni1.^2 .* exp(1i * (phase1 - p0cs)), 4));
combined2 = angle(sum(magni2.^2 .* exp(1i * (phase2 - p0cs)), 4));
end


% low pass filter the estimate of coils phase offset
function lp = lp_filter(magni, phase, spatial_res)

cmplx = magni.*exp(1i*phase);
rimg = real(cmplx);
iimg = imag(cmplx);

rSmooth = rimg;
fprintf('\nSmoothing coil: ')
for i = 1:size(rimg, 4)
    fprintf('%d ', i);
    rSmooth(:,:,:,i) = imgaussfilt3(rimg(:,:,:,i), 1.5./spatial_res);
end

iSmooth = iimg;
fprintf('\nSmoothing coil: ')
for i = 1:size(iimg, 4)
    fprintf('%d ', i);
    iSmooth(:,:,:,i) = imgaussfilt3(iimg(:,:,:,i), 1.5./spatial_res);
end
fprintf('\n')

lp = angle(rSmooth + 1i * iSmooth);

end

function eroded = erode_mask(mask, R)
[xx,yy,zz] = ndgrid(-R:R);
nhood = sqrt(xx.^2 + yy.^2 + zz.^2) <= R;
eroded = imerode(mask,nhood);
end
