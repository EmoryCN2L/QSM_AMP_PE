function hdr = niftiHdr();
%niftiHdr create an empty nifti header
%hdr.hk, hdr.dime, hdr.hist
%Key elements to fill in
%   hdr.dime.dim
%   hdr.dime.pixdim
%   hdr.dime.datatype
%   hdr.dime.scl_inter, hdr.dime.scl_slope

    %  Original header structures
	%  struct dsr
	%       { 
	%       struct header_key hk;            /*   0 +  40       */
	%       struct image_dimension dime;     /*  40 + 108       */
	%       struct data_history hist;        /* 148 + 200       */
	%       };                               /* total= 348 bytes*/

    hdr.hk   = header_key;
    hdr.dime = image_dimension;
    hdr.hist = data_history;

    %  For Analyze data format
    %
    if ~strcmp(hdr.hist.magic, 'n+1') & ~strcmp(hdr.hist.magic, 'ni1')
        dsr.hist.qform_code = 0;
        dsr.hist.sform_code = 0;
    end

    return					% read_header


%---------------------------------------------------------------------
function [ hk ] = header_key();

    
	%  Original header structures	
	%  struct header_key                     /* header key      */ 
	%       {                                /* off + size      */
	%       int sizeof_hdr                   /*  0 +  4         */
	%       char data_type[10];              /*  4 + 10         */
	%       char db_name[18];                /* 14 + 18         */
	%       int extents;                     /* 32 +  4         */
	%       short int session_error;         /* 36 +  2         */
	%       char regular;                    /* 38 +  1         */
	%       char dim_info;   % char hkey_un0;        /* 39 +  1 */
	%       };                               /* total=40 bytes  */
	%
	% int sizeof_header   Should be 348.
	% char regular        Must be 'r' to indicate that all images and 
	%                     volumes are the same size. 

    v6 = version;
    if str2num(v6(1))<6
       directchar = '*char';
    else
       directchar = 'uchar=>char';
    end

    hk.sizeof_hdr    = 348;% should be 348!
    hk.data_type     = 'float32';
    hk.db_name       = '';
    hk.extents       = 0;
    hk.session_error = 0;
    hk.regular       = 'r';
    hk.dim_info      = 0;
    
    return					% header_key


%---------------------------------------------------------------------
function [ dime ] = image_dimension()

	%  Original header structures    
	%  struct image_dimension
	%       {                                /* off + size      */
	%       short int dim[8];                /* 0 + 16          */
        %       /*
        %           dim[0]      Number of dimensions in database; usually 4. 
        %           dim[1]      Image X dimension;  number of *pixels* in an image row. 
        %           dim[2]      Image Y dimension;  number of *pixel rows* in slice. 
        %           dim[3]      Volume Z dimension; number of *slices* in a volume. 
        %           dim[4]      Time points; number of volumes in database
        %       */
	%       float intent_p1;   % char vox_units[4];   /* 16 + 4       */
	%       float intent_p2;   % char cal_units[8];   /* 20 + 4       */
	%       float intent_p3;   % char cal_units[8];   /* 24 + 4       */
	%       short int intent_code;   % short int unused1;   /* 28 + 2 */
	%       short int datatype;              /* 30 + 2          */
	%       short int bitpix;                /* 32 + 2          */
	%       short int slice_start;   % short int dim_un0;   /* 34 + 2 */
	%       float pixdim[8];                 /* 36 + 32         */
	%	/*
	%		pixdim[] specifies the voxel dimensions:
	%		pixdim[1] - voxel width, mm
	%		pixdim[2] - voxel height, mm
	%		pixdim[3] - slice thickness, mm
	%		pixdim[4] - volume timing, in msec
	%					..etc
	%	*/
	%       float vox_offset;                /* 68 + 4          */
	%       float scl_slope;   % float roi_scale;     /* 72 + 4 */
	%       float scl_inter;   % float funused1;      /* 76 + 4 */
	%       short slice_end;   % float funused2;      /* 80 + 2 */
	%       char slice_code;   % float funused2;      /* 82 + 1 */
	%       char xyzt_units;   % float funused2;      /* 83 + 1 */
	%       float cal_max;                   /* 84 + 4          */
	%       float cal_min;                   /* 88 + 4          */
	%       float slice_duration;   % int compressed; /* 92 + 4 */
	%       float toffset;   % int verified;          /* 96 + 4 */
	%       int glmax;                       /* 100 + 4         */
	%       int glmin;                       /* 104 + 4         */
	%       };                               /* total=108 bytes */
	
    dime.dim        = [4,ones(1,7)];
    dime.intent_p1  = 0;
    dime.intent_p2  = 0;
    dime.intent_p3  = 0;
    dime.intent_code = 0;
    dime.datatype   = niftiType('float32');
    dime.bitpix     = 32;
    dime.slice_start = 0;
    dime.pixdim     = [-1,ones(1,7)];
    dime.vox_offset = 352;
    dime.scl_slope  = 1;
    dime.scl_inter  = 0;
    dime.slice_end  = 0;
    dime.slice_code = 0;
    dime.xyzt_units = 10;
    dime.cal_max    = 0;
    dime.cal_min    = 0;
    dime.slice_duration = 0;
    dime.toffset    = 0;
    dime.glmax      = 0;
    dime.glmin      = 0;
        
    return					% image_dimension


%---------------------------------------------------------------------
function [ hist ] = data_history()
        
	%  Original header structures
	%  struct data_history       
	%       {                                /* off + size      */
	%       char descrip[80];                /* 0 + 80          */
	%       char aux_file[24];               /* 80 + 24         */
	%       short int qform_code;            /* 104 + 2         */
	%       short int sform_code;            /* 106 + 2         */
	%       float quatern_b;                 /* 108 + 4         */
	%       float quatern_c;                 /* 112 + 4         */
	%       float quatern_d;                 /* 116 + 4         */
	%       float qoffset_x;                 /* 120 + 4         */
	%       float qoffset_y;                 /* 124 + 4         */
	%       float qoffset_z;                 /* 128 + 4         */
	%       float srow_x[4];                 /* 132 + 16        */
	%       float srow_y[4];                 /* 148 + 16        */
	%       float srow_z[4];                 /* 164 + 16        */
	%       char intent_name[16];            /* 180 + 16        */
	%       char magic[4];   % int smin;     /* 196 + 4         */
	%       };                               /* total=200 bytes */

    v6 = version;
    if str2num(v6(1))<6
       directchar = '*char';
    else
       directchar = 'uchar=>char';
    end

    hist.descrip     = '';
    hist.aux_file    = '';
    hist.qform_code  = 0;
    hist.sform_code  = 0;
    hist.quatern_b   = 0;
    hist.quatern_c   = 0;
    hist.quatern_d   = 0;
    hist.qoffset_x   = 0;
    hist.qoffset_y   = 0;
    hist.qoffset_z   = 0;
    hist.srow_x      = [1,0,0,0];
    hist.srow_y      = [0,1,0,0];
    hist.srow_z      = [0,0,1,0];
    hist.intent_name = '';
    hist.magic       = 'n+1';

  
    hist.originator  = zeros(1,5);
    
    return					% data_history

