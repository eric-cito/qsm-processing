function hdr = read_raw_hdr_26x_Jingwen(pfilename)
% read_raw_hdr info 
% ---------------------
% Usage: >> hdr = read_raw_hdr('pfile_name')
% Description: Reads the pfile header values necessary for reconstructing
% the pfile image and creating a .idf file%
%
% Modified from read_raw_hdr_new to work with svk_gepfile_reader.dev (dv26)
% Eason Chen, 09/01/2018
%

disp(['read raw header pfile name: ' pfilename]);

if (nargin < 1)
    error('''pfile_info.m'' was not passed a P-file name');
end

field_defs = ...
    {'hdr_rec.rdb_hdr_rdbm_rev =',         '%f',  inf, 'rdb_hdr_rdbm_rev'    , []
    'hdr_rec.rdb_hdr_da_xres =',           '%d',  inf, 'rdb_hdr_da_xres',  []
    'hdr_rec.rdb_hdr_da_yres =',           '%d',  inf, 'rdb_hdr_da_yres',  []
    'hdr_rec.rdb_hdr_rc_xres =',           '%d',  inf, 'rdb_hdr_rc_xres',  []
    'hdr_rec.rdb_hdr_rc_yres =',           '%d',  inf, 'rdb_hdr_rc_yres',  []
    'hdr_rec.rdb_hdr_frame_size =',        '%d',  inf, 'rdb_hdr_frame_size',  []
    'hdr_rec.rdb_hdr_nframes =',           '%d',  inf, 'rdb_hdr_nframes',  []
    'hdr_rec.rdb_hdr_im_size =',           '%d',  inf, 'rdb_hdr_im_size',  []
    'hdr_rec.rdb_hdr_te =',                '%d',  inf, 'rdb_hdr_te',       []
    'hdr_rec.rdb_hdr_te2 =',               '%d',  inf, 'rdb_hdr_te2',      []
    'hdr_rec.rdb_hdr_scalei =',            '%d',  inf, 'rdb_hdr_scalei',  []
    'hdr_rec.rdb_hdr_point_size =',        '%d',  inf, 'rdb_hdr_point_size',  []
    'hdr_rec.rdb_hdr_raw_pass_size =',     '%d',  inf, 'rdb_hdr_raw_pass_size',  []
    'hdr_rec.rdb_hdr_phase_scale =',       '%d',  inf, 'rdb_hdr_phase_scale',  []
    'hdr_rec.rdb_hdr_nslices =',           '%d',  inf, 'rdb_hdr_nslices',  []
    'hdr_rec.rdb_hdr_slblank =',           '%d',  inf, 'rdb_hdr_slblank',  []
    'hdr_rec.rdb_hdr_ileaves =',           '%d',  inf, 'rdb_hdr_ileaves',  []
    'hdr_rec.rdb_hdr_navs =',              '%d',  inf, 'rdb_hdr_navs',  [] 
    'hdr_rec.rdb_hdr_nechoes =',           '%d',  inf, 'rdb_hdr_nechoes',  []
    'hdr_rec.rdb_hdr_npasses =',           '%d',  inf, 'rdb_hdr_npasses',  []
    'hdr_rec.rdb_hdr_fov =',               '%d',  inf, 'rdb_hdr_fov',  []
    'hdr_rec.rdb_hdr_scancent =',          '%f',  inf, 'rdb_hdr_scancent',  []       
    'hdr_rec.rdb_hdr_data_collect_type =', '%d',  inf, 'rdb_hdr_data_collect_type',  []
    'hdr_rec.rdb_hdr_recon_ctrl =',        '%d',  inf, 'rdb_hdr_recon_ctrl',  []
    'hdr_rec.rdb_hdr_exec_ctrl =',         '%d',  inf, 'rdb_hdr_exec_ctrl',  []
    'hdr_rec.rdb_hdr_dacq_ctrl =',         '%d',  inf, 'rdb_hdr_dacq_ctrl',  []
    'hdr_rec.rdb_hdr_scan_date =',         '%s',  inf, 'rdb_hdr_scan_date',  []
    'hdr_rec.rdb_hdr_dab[0].start_rcv =',  '%d',  inf, 'rdb_hdr_start_rcv',  []
    'hdr_rec.rdb_hdr_dab[0].stop_rcv =',   '%d',  inf, 'rdb_hdr_stop_rcv',  []
    'hdr_rec.rdb_hdr_ovl =',               '%d',  inf, 'rdb_hdr_ovl',  []
    'exam.ex_no =',                        '%s',  inf, 'ex_no',  []
    'exam.patid =',                        '%s',  inf, 'patid',  []
    'exam.patname =',                      '%s',  inf, 'patname',  []
    'series.se_desc =',                    '%c',  inf, 'se_desc',  [] 
    'series.se_no =',                      '%d',  inf, 'se_no',  [] 
    'series.se_exno =',                    '%d',  inf, 'se_exno',  []
    'series.entry ='                       '%d',  inf, 'entry',  []
    'series.position =',                   '%d',  inf, 'position',  []
    'series.start_loc =',                  '%f',  inf, 'start_loc',  []
    'series.end_loc =',                    '%f',  inf, 'end_loc',  []
    'image.cname =',                       '%c',  inf, 'cname', []
    'image.imode =',                       '%d',  inf, 'imode',  []
    'image.imatrix_X =',                    '%d',  inf, 'imatrix_X',  []
    'image.imatrix_Y =',                    '%d',  inf, 'imatrix_Y',  []
    'image.freq_dir =',                    '%d',  inf, 'freq_dir',  []
    'image.fphase =',                      '%d',  inf, 'fphase',  []
    'image.plane =',                       '%d',  inf, 'plane',  []
    'image.psdname =',                     '%c',  inf, 'psdname',  []
    'image.dim_X =',                       '%d',  inf, 'dim_X',  []
    'image.dim_Y =',                       '%d',  inf, 'dim_Y',  []
    'image.slthick =',                     '%f',  inf, 'slthick',  []
    'image.scanspacing =',                 '%f',  inf, 'scanspacing',  []
    'image.numslabs =',                    '%f',  inf, 'numslabs',  []
    'image.locsperslab =',                 '%f',  inf, 'locsperslab',  []
    'image.ctr_R =',                       '%f',  inf, 'ctr_R',  []
    'image.ctr_A =',                       '%f',  inf, 'ctr_A',  []
    'image.ctr_S =',                       '%f',  inf, 'ctr_S',  []
    'image.slthick =',                     '%f',  inf, 'slthick',  []
    'image.tlhc_R =',                      '%f',  inf, 'tlhc_R',  []
    'image.tlhc_A =',                      '%f',  inf, 'tlhc_A',  []
    'image.tlhc_S =',                      '%f',  inf, 'tlhc_S',  []
    'image.trhc_R =',                      '%f',  inf, 'trhc_R',  []
    'image.trhc_A =',                      '%f',  inf, 'trhc_A',  []
    'image.trhc_S =',                      '%f',  inf, 'trhc_S',  []
    'image.brhc_R =',                      '%f',  inf, 'brhc_R',  []
    'image.brhc_A =',                      '%f',  inf, 'brhc_A',  []
    'image.brhc_S =',                      '%f',  inf, 'brhc_S',  []
    'image.norm_R =',                      '%f',  inf, 'norm_R',  []
    'image.norm_A =',                      '%f',  inf, 'norm_A',  []
    'image.norm_S =',                      '%f',  inf, 'norm_S',  []
    'image.dfov =',                        '%f',  inf, 'dfov',    []        % X
    'image.dfov_rect =',                   '%f',  inf, 'dfov_rect', []      % Y      
};



fid = fopen(pfilename, 'r', 'l');
version = fread(fid, 1, 'float');
fclose(fid);

% old pfile
[status entire_file] = system(sprintf('rdump.dev -Ap %s', pfilename));

disp(['read raw header version: ' num2str(version)]);
if (version >= 11)
    field_defs_14x={'hdr_rec.rdb_hdr_off_data =',          '%d',  inf, 'rdb_hdr_off_data',  [] 
                    'hdr_rec.rdb_hdr_off_data_acq_tab =',  '%d',  inf, 'rdb_hdr_off_data_acq_tab', [] };
    field_defs=cat(1, field_defs, field_defs_14x);
end


% dv26
if (strfind(entire_file, 'Unsupported version of raw header'))
    [status entire_file] = system(sprintf('svk_gepfile_reader.dev -i %s --print_header', pfilename));
    entire_file = regexprep(entire_file, ' +',' ');
    field_defs = ...
        {'rhr.rh_rdbm_rev =',         '%f',  inf, 'rdb_hdr_rdbm_rev'    , []
        'rhr.rh_da_xres =',           '%d',  inf, 'rdb_hdr_da_xres',  []
        'rhr.rh_da_yres =',           '%d',  inf, 'rdb_hdr_da_yres',  []
        'rhr.rh_rc_xres =',           '%d',  inf, 'rdb_hdr_rc_xres',  []
        'rhr.rh_rc_yres =',           '%d',  inf, 'rdb_hdr_rc_yres',  []
        'rhr.rh_frame_size =',        '%d',  inf, 'rdb_hdr_frame_size',  []
        'rhr.rh_nframes =',           '%d',  inf, 'rdb_hdr_nframes',  []
        'rhr.rh_im_size =',           '%d',  inf, 'rdb_hdr_im_size',  []
        'rhr.rh_te =',                '%d',  inf, 'rdb_hdr_te',       []
        'rhr.rh_te2 =',               '%d',  inf, 'rdb_hdr_te2',      []
        'rhr.rh_scalei =',            '%d',  inf, 'rdb_hdr_scalei',  []
        'rhr.rh_point_size =',        '%d',  inf, 'rdb_hdr_point_size',  []
        'rhr.rh_raw_pass_size =',     '%d',  inf, 'rdb_hdr_raw_pass_size',  []
        'rhr.rh_phase_scale =',       '%d',  inf, 'rdb_hdr_phase_scale',  []
        'rhr.rh_nslices =',           '%d',  inf, 'rdb_hdr_nslices',  []
        'rhr.rh_slblank =',           '%d',  inf, 'rdb_hdr_slblank',  []
        'rhr.rh_ileaves =',           '%d',  inf, 'rdb_hdr_ileaves',  []
        'rhr.rh_navs =',              '%d',  inf, 'rdb_hdr_navs',  [] 
        'rhr.rh_nechoes =',           '%d',  inf, 'rdb_hdr_nechoes',  []
        'rhr.rh_npasses =',           '%d',  inf, 'rdb_hdr_npasses',  []
        'rhr.rh_fov =',               '%d',  inf, 'rdb_hdr_fov',  []
        'rhr.rh_scancent =',          '%f',  inf, 'rdb_hdr_scancent',  []       
        'rhr.rh_data_collect_type =', '%d',  inf, 'rdb_hdr_data_collect_type',  []
        'rhr.rh_recon_ctrl =',        '%d',  inf, 'rdb_hdr_recon_ctrl',  []
        'rhr.rh_exec_ctrl =',         '%d',  inf, 'rdb_hdr_exec_ctrl',  []
        'rhr.rh_dacq_ctrl =',         '%d',  inf, 'rdb_hdr_dacq_ctrl',  []
        'rhr.rh_scan_date =',         '%s',  inf, 'rdb_hdr_scan_date',  []
        'rhr.rh_dab[0].start_rcv =',  '%d',  inf, 'rdb_hdr_start_rcv',  []
        'rhr.rh_dab[0].stop_rcv =',   '%d',  inf, 'rdb_hdr_stop_rcv',  []
        'rhr.rh_ovl =',               '%d',  inf, 'rdb_hdr_ovl',  []
        'rhe.ex_no =',                        '%s',  inf, 'ex_no',  []
        'rhe.patid =',                        '%s',  inf, 'patid',  []
        'rhe.patname =',                      '%s',  inf, 'patname',  []
        'rhs.se_desc =',                    '%c',  inf, 'se_desc',  [] 
        'rhs.se_no =',                      '%d',  inf, 'se_no',  [] 
        'rhs.se_exno =',                    '%d',  inf, 'se_exno',  []
        'rhs.entry ='                       '%d',  inf, 'entry',  []
        'rhs.position =',                   '%d',  inf, 'position',  []
        'rhs.start_loc =',                  '%f',  inf, 'start_loc',  []
        'rhs.end_loc =',                    '%f',  inf, 'end_loc',  []
        'rhi.cname =',                       '%c',  inf, 'cname', []
        'rhi.imode =',                       '%d',  inf, 'imode',  []
        'rhi.imatrix_X =',                    '%d',  inf, 'imatrix_X',  []
        'rhi.imatrix_Y =',                    '%d',  inf, 'imatrix_Y',  []
        'rhi.freq_dir =',                    '%d',  inf, 'freq_dir',  []
        'rhi.fphase =',                      '%d',  inf, 'fphase',  []
        'rhi.plane =',                       '%d',  inf, 'plane',  []
        'rhi.psdname =',                     '%c',  inf, 'psdname',  []
        'rhi.dim_X =',                       '%d',  inf, 'dim_X',  []
        'rhi.dim_Y =',                       '%d',  inf, 'dim_Y',  []
        'rhi.slthick =',                     '%f',  inf, 'slthick',  []
        'rhi.scanspacing =',                 '%f',  inf, 'scanspacing',  []
        'rhi.numslabs =',                    '%f',  inf, 'numslabs',  []
        'rhi.locsperslab =',                 '%f',  inf, 'locsperslab',  []
        'rhi.ctr_R =',                       '%f',  inf, 'ctr_R',  []
        'rhi.ctr_A =',                       '%f',  inf, 'ctr_A',  []
        'rhi.ctr_S =',                       '%f',  inf, 'ctr_S',  []
        'rhi.slthick =',                     '%f',  inf, 'slthick',  []
        'rhi.tlhc_R =',                      '%f',  inf, 'tlhc_R',  []
        'rhi.tlhc_A =',                      '%f',  inf, 'tlhc_A',  []
        'rhi.tlhc_S =',                      '%f',  inf, 'tlhc_S',  []
        'rhi.trhc_R =',                      '%f',  inf, 'trhc_R',  []
        'rhi.trhc_A =',                      '%f',  inf, 'trhc_A',  []
        'rhi.trhc_S =',                      '%f',  inf, 'trhc_S',  []
        'rhi.brhc_R =',                      '%f',  inf, 'brhc_R',  []
        'rhi.brhc_A =',                      '%f',  inf, 'brhc_A',  []
        'rhi.brhc_S =',                      '%f',  inf, 'brhc_S',  []
        'rhi.norm_R =',                      '%f',  inf, 'norm_R',  []
        'rhi.norm_A =',                      '%f',  inf, 'norm_A',  []
        'rhi.norm_S =',                      '%f',  inf, 'norm_S',  []
        'rhi.dfov =',                        '%f',  inf, 'dfov',    []        % X
        'rhi.dfov_rect =',                   '%f',  inf, 'dfov_rect', []      % Y      
    };

    if (version >= 11)
        field_defs_14x={'rhr.rdb_hdr_off_data =',          '%d',  inf, 'rdb_hdr_off_data',  []};
%                         'rhr.rdb_hdr_off_data_acq_tab =',  '%d',  inf, 'rdb_hdr_off_data_acq_tab', [] };
        field_defs=cat(1, field_defs, field_defs_14x);
    end
end

if (numel(entire_file) < 1000)
    disp(['PFILE (name, numelements): ' pfilename, ' ', num2str(numel(entire_file))]);
end

number_of_fields = size(field_defs,1);

location_end = 1;

hdr = [];

% =======================================================================
% Begin reading the actual data
% =======================================================================

for (fieldnum = 1:number_of_fields)
    
%     if fieldnum == 72
%         disp('1');
%     end

    % Find the field (assume line length < 500)
    location_start = strfind(entire_file, field_defs{fieldnum,1});
    location_end = strfind(entire_file(location_start:min(location_start + 500, length(entire_file))), sprintf('\n'));
    field_name = field_defs{fieldnum,4};
    
%     fprintf('Test - Jingwen: %d-%d, %d\n', location_start, location_end, length(entire_file));
%     fprintf('%d-%d, %s\n', location_start, location_end(1), field_name);
    
    field_value = sscanf(entire_file(location_start(1)+length(field_defs{fieldnum,1}):location_start(1) + location_end(1)-1), ...
        field_defs{fieldnum,2},                   ...
        field_defs{fieldnum,3});
    
    % If it's a string, remove white space before and after the desired word
    if (isstr(field_value))
        delim = strfind(field_value, '#');
        if (length(delim))
            field_value = field_value(1:delim(1)-1);
        end
        field_value = strtrim(field_value);
    end
    
    % Assign the data into the structure
    
    array_index = field_defs{fieldnum,5};
    
    if (isempty(array_index))
        hdr = setfield(hdr, field_name, field_value);
    else
        if (isfield(pfile, field_name))
            temp_value = getfield(hdr, field_name);
        else
            temp_value = [];
        end
        if (isstr(field_value))
            temp_value{array_index} = field_value;
        else
            temp_value(array_index) = field_value;
        end
        hdr = setfield(hdr, field_name, temp_value);
    end
  
end

% =======================================================================
% Added for reading data acq_tab section of the header .... necessary for
% combining data from multiple Pfiles
% =======================================================================

clear fieldnum field_value field_name

new_field_defs=[];
nslice=hdr.rdb_hdr_nslices;
im_mode=hdr.imode;

if (im_mode == 1)

    for (n=1:nslice)
        passnum_strname=strcat('acq_tab[', num2str(n-1), '].pass_number =');
        slicepass_strname=strcat('acq_tab[', num2str(n-1), '].slice_in_pass =');
        field_defs_acq_tab={passnum_strname,          '%d',  inf, 'pass_number',  [n]
            slicepass_strname,        '%d',  inf, 'slice_in_pass',  [n]   };
        new_field_defs=cat(1, new_field_defs, field_defs_acq_tab);
    end
    number_of_fields = size(new_field_defs,1);


    for (fieldnum = 1:number_of_fields)

        % Find the field (assume line length < 500)
        location_start = strfind(entire_file, new_field_defs{fieldnum,1});
        location_end = strfind(entire_file(location_start:location_start + 500), sprintf('\n'));
        field_name  = new_field_defs{fieldnum,4};
        field_value = sscanf(entire_file(location_start(1)+length(new_field_defs{fieldnum,1}):location_start(1) + location_end(1)-1), ...
            new_field_defs{fieldnum,2},                   ...
            new_field_defs{fieldnum,3});

        % If it's a string, remove white space before and after the desired word

        if (isstr(field_value))
            delim = strfind(field_value, '#');
            if (length(delim))
                field_value = field_value(1:delim(1)-1);
            end
            field_value = strtrim(field_value);
        end

        % Assign the data into the structure

        array_index = new_field_defs{fieldnum,5};

        if (isempty(array_index))
            hdr = setfield(hdr, field_name, field_value);
        else
            if (isfield(hdr, field_name))
                temp_value = getfield(hdr, field_name);
            else
                temp_value = [];
            end
            if (isstr(field_value))
                temp_value{array_index} = field_value;
            else
                temp_value(array_index) = field_value;
            end

            hdr = setfield(hdr, field_name, temp_value);
        end

    end

end


