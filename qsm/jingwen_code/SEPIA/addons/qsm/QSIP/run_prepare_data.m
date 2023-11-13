function [atlas, noisemask, est0] = run_prepare_data(pn_pre, target_dim)
fn_est0 = 'susc_est_ppm_0';
fn_atlas = 'ROR_atlas';
fn_noisemask = 'ROR_noise_mask';

[img,dims,scales,bpp,endian] = read_avw([pn_pre fn_atlas]);
atlas = zeros([target_dim size(img, 4)]);
for i = 1:size(img, 4)
    atlas(:, :, :, i) = imresize3(img(:, :, :, i), target_dim);
end

[img,dims,scales,bpp,endian] = read_avw([pn_pre fn_noisemask]);
noisemask = zeros([target_dim size(img, 4)]);
for i = 1:size(img, 4)
    noisemask(:, :, :, i) = imresize3(img(:, :, :, i), target_dim);
end

[img,dims,scales,bpp,endian] = read_avw([pn_pre fn_est0]);
est0 = zeros([target_dim size(img, 4)]);
for i = 1:size(img, 4)
    est0(:, :, :, i) = imresize3(img(:, :, :, i), target_dim);
end