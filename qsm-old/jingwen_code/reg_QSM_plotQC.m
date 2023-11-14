function [] = reg_QSM_plotQC(QSMname, ROIname, targetSize, flagflip)

if nargin < 4
    flagflip = 1;
end

nii = load_untouch_nii(QSMname);
QSM = nii.img;
currSize = size(QSM);
if currSize(1) > targetSize(1)
    targetSize(1:2) = currSize(1);
end
if currSize(2) > targetSize(2)
    targetSize(1:2) = currSize(2);
end
if currSize(3) > targetSize(3)
    targetSize(3) = currSize(3);
end
padSize = (targetSize - currSize)/2;

midSize = floor(targetSize/2);

nii = load_untouch_nii(ROIname);
QSMseg = nii.img;

if flagflip
    QSM_Ax = padarray(flip(rot90(QSM(:,:,midSize(3)),1)),padSize([2 1]));
    QSM_Cor = padarray(flip(rot90(squeeze(QSM(midSize(1),:,:)),1)),padSize([3 2]));
    QSM_Sag = padarray(flip(rot90(squeeze(QSM(:,midSize(2),:)),1)),padSize([3 1]));
    QSMseg_Ax = padarray(flip(rot90(QSMseg(:,:,midSize(3)),1)),padSize([2 1]));
    QSMseg_Cor = padarray(flip(rot90(squeeze(QSMseg(midSize(1),:,:)),1)),padSize([3 2]));
    QSMseg_Sag = padarray(flip(rot90(squeeze(QSMseg(:,midSize(2),:)),1)),padSize([3 1]));
else
    QSM_Ax = padarray((rot90(QSM(:,:,midSize(3)),1)),padSize([2 1]));
    QSM_Cor = padarray((rot90(squeeze(QSM(midSize(1),:,:)),1)),padSize([3 2]));
    QSM_Sag = padarray((rot90(squeeze(QSM(:,midSize(2),:)),1)),padSize([3 1]));
    QSMseg_Ax = padarray((rot90(QSMseg(:,:,midSize(3)),1)),padSize([2 1]));
    QSMseg_Cor = padarray((rot90(squeeze(QSMseg(midSize(1),:,:)),1)),padSize([3 2]));
    QSMseg_Sag = padarray((rot90(squeeze(QSMseg(:,midSize(2),:)),1)),padSize([3 1]));
end

X = [QSM_Ax; QSM_Cor; QSM_Sag];
Y = [QSMseg_Ax; QSMseg_Cor; QSMseg_Sag];

figure('position',[100 100 500 500]); 
ax1 = axes; imagesc(X); colormap(ax1,'gray'); % [0 0.3]
axis tight; axis off; axis equal;
ax2 = axes; imagesc(Y, 'alphadata', (Y>0)*0.3); colormap(ax2,'jet');
caxis(ax2,[min(nonzeros(Y)) max(nonzeros(Y))]);
ax2.Visible = 'off';
axis tight; axis off; axis equal;
linkprop([ax1 ax2],'Position');

end