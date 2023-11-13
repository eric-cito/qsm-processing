
%% load data
FileStruct = load_untouch_nii('PDa391_qsm_e1.nii');
Mag1 = FileStruct.img;
FileStruct = load_untouch_nii('PDa391_qsm_e1_ph.nii');
Ph1 = FileStruct.img;


%% scale and attempt to fix phase chop
% minPh = min(Ph1(:)); maxPh = max(Ph1(:));
% rangePh = maxPh - minPh;
%tmp_ph = double(Ph1-minPh)/double(rangePh)*2*pi; %changing phase range
tmp_ph = double(Ph1)/4096*pi;
image; imagesc(tmp_ph(:,:,72),[-pi pi]); colorbar
compl = double(Mag1).*exp(1i*tmp_ph); %combining mag and phase, complex
compl_re = ifft(ifftshift(fft(compl,[],3),3),[],3);
compl_k = fftdim(compl_re,[1 2 3]); %fft to create k space

figure;
subplot(1,3,1); imagesc(abs(compl_k(:,:,72)).^0.2); %k-space image
subplot(1,3,2); imagesc(abs(squeeze(compl_k(:,256,:))).^0.2); %k-space
subplot(1,3,3); imagesc(abs(squeeze(compl_k(256,:,:))).^0.2); %k-space

%% apply  hamming in attempt to further fix phase chop

w1 = hamming(size(compl,1)); %hamming filter
w2 = hamming(size(compl,3)); %hamming filter
w3 = repmat(reshape(w1*w2',1,size(compl,1),[]),[size(compl,1) 1 1]);
compl_k = compl_k.*w3;

figure;
subplot(1,3,1); imagesc(abs(compl_k(:,:,72)).^0.2); %k-space image
subplot(1,3,2); imagesc(abs(squeeze(compl_k(:,256,:))).^0.2); %k-space
subplot(1,3,3); imagesc(abs(squeeze(compl_k(256,:,:))).^0.2); %k-space

compl_re = ifftdim(compl_k,[1 2 3]);

figure; 
subplot(1,2,1); 
imagesc(angle(squeeze(compl_re(256,:,:))));
subplot(1,2,2); 
imagesc(abs(squeeze(compl_re(256,:,:))));

%% save image

FileStruct.img=angle(compl_re);
FileStruct.hdr.dime.cal_max= pi;
FileStruct.hdr.dime.cal_min=-pi;
FileStruct.hdr.dime.datatype=64;
FileStruct.hdr.dime.bitpix=16;
save_untouch_nii(FileStruct, 'PDa391_fix_e1_ph.nii');

