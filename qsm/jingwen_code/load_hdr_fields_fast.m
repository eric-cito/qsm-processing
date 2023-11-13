%  Usage: load_hdr_fields
%
%  loads header info obtained from pfile header
%  created from the structure 'hdr' output from
%  read_pfile_dev and populates variables used 
%  in most recon code
%
%  For example:
%    >> hdr=read_pfile_dev('P123456.7');
%    >> load_hdr_fields  
%
%  Janine Lupo
%  April 2007
%  University of California, San Francisco
%

%frequency
xdim=hdr.rdb_hdr_da_xres;
xres=hdr.rdb_hdr_frame_size;
im_xres=hdr.rdb_hdr_rc_xres;
scan_freq=hdr.dim_X;

%phase
ydim=hdr.rdb_hdr_da_yres;
yres=hdr.rdb_hdr_nframes;
im_yres=hdr.rdb_hdr_rc_yres;
scan_phase=hdr.dim_Y;

pfov=hdr.rdb_hdr_phase_scale; %% pfov= phase field of view
scan_phase=ceil(pfov*scan_phase);   
if rem(scan_phase,2)~=0
    scan_phase=scan_phase+1;
end

if (hdr.dfov == hdr.dfov_rect) || (hdr.dfov_rect == 0)
    fov=hdr.rdb_hdr_fov;    % square FOV
else
    fov = [hdr.dfov_rect hdr.dfov];
end

im_size=hdr.rdb_hdr_im_size;

im_mode=hdr.imode;  % 1 for 2D, 2 for 3D
freq_dir=hdr.freq_dir; % 1 for R/L, 2 for A/P
%bestky=hdr.rdb_hdr_pcbestky;

% determine rhrcctrl, rhexecctrl, and rhtype bitmasks 
% and set appropriate flags

%scale_factor=hdr.rdb_hdr_scalei;
rhrcctrl=hdr.rdb_hdr_recon_ctrl;
rhexecctrl=hdr.rdb_hdr_exec_ctrl;
rhdacqctrl=hdr.rdb_hdr_dacq_ctrl;
rhtype=hdr.rdb_hdr_data_collect_type;

if bitget(rhtype, 1) == 1
    chop_flag=1;
else
    chop_flag=0;
end

if bitget(rhtype, 9) == 1   
    npw=1;  %no phase wrap on
else
    npw=0;
end

if bitget(rhdacqctrl, 1) == 1
    rawdata=1;
else
    rawdata=0;
end



if bitget(rhrcctrl, 15) == 1   
    ifft3d_flag=1;
else
    ifft3d_flag=0;
end

%slice (or 2nd phase)

nslice=hdr.rdb_hdr_nslices;
nblank=hdr.rdb_hdr_slblank;
npass=hdr.rdb_hdr_npasses;
sl_pass=hdr.rdb_hdr_nslices/npass;
%acq_tab=hdr.acq_tab;

if im_mode==1
    first_slice=1;
    last_slice=nslice;
elseif im_mode==2 || im_mode==9
    nslice=sl_pass-2*nblank;
    first_slice = nblank+1;
    last_slice = sl_pass-nblank;
else
    nslice=1;
    im_mode==1;
end

%added fields for multi-slab
novl=hdr.rdb_hdr_ovl;
nlocslab=hdr.locsperslab;
nslab=hdr.numslabs;

%extra dimensions
nex=hdr.rdb_hdr_navs;
necho=hdr.rdb_hdr_nechoes;

if hdr.fphase>=1
    nphase=hdr.fphase;
else
    nphase=1;
end


%view (y) interleave for multishot data
ileaves=hdr.rdb_hdr_ileaves;
if ileaves ~= 0
    vpsht=yres/ileaves;
end


%determining precision
pt_size=hdr.rdb_hdr_point_size;  % 2 for regular, 4 for EDR 
if pt_size==2
    datatype='int16';
    bitsize=16;
elseif pt_size==4
    datatype='int32';
    bitsize=32;
else
    sprintf('Invalid pointsize in header.  Using int16 datatype.');
    datatype='int16';
end

%determine header size
rdbm_version=hdr.rdb_hdr_rdbm_rev;
if rdbm_version==9       % 11x 
    hdr_size=61464;
elseif rdbm_version==11  % 12x
    hdr_size=66072;  
else                     % 14x
    hdr_size=hdr.rdb_hdr_off_data;  
end

%determining # of coils from # of values in receiver noise matrix  or
%filesize 
% if (hdr.rdb_hdr_rdbm_rev<14)
%     ncoil=length(find(hdr.rdb_hdr_rec_noise_std~=1));
%     noise_std=hdr.rdb_hdr_rec_noise_std(find(hdr.rdb_hdr_rec_noise_std~=1));
% else
%    ncoil=hdr.rdb_hdr_raw_pass_size./(xdim*(ydim)*sl_pass*npass*necho*pt_size*2);
%end

ncoil = hdr.rdb_hdr_stop_rcv-hdr.rdb_hdr_start_rcv +1;

%check to make sure rawpass size equals the correct # of bytes
rawpass_size=hdr.rdb_hdr_raw_pass_size;
temp_rawpass=xdim*(ydim)*sl_pass*necho*ncoil*pt_size*2;

if rawpass_size~=temp_rawpass
    sprintf('Warning: File size does not match dimensions in header');
end

raw_vol=xdim*(ydim)*sl_pass*ncoil*2;
coil_vol=xdim*(ydim)*sl_pass*2;
slice_vol=xdim*(ydim)*2;

if ncoil < 1
    ncoil=1;
end

%to force S/I coordinates for 2D acquisitions only
if (hdr.end_loc <=  hdr.start_loc)    
    SI_flag=0;  %correct
else
    SI_flag=1;  %need to flip data in z
end
   

