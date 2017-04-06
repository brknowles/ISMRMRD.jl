module ISMRMRD

__precompile__(true)

import Base: start, next, done, length, size, eltype;

# lots to export
export 	ISMRMRD_dataset,
	ISMRMRD_acquisition,
	ISMRMRD_acquisitionheader,
	ISMRMRD_EncodingCounters,
	read_xml_header,
	write_xml_header,
	append_dataset,
	read_acquisition,
	get_number_of_acquisitions,
	get_number_of_data_samples,
	get_number_of_trajectory_samples,
	get_data,
	get_trajectory,
	get_version,
	get_flags,
	get_meas_uid,
	get_scan_counter,
	get_acquisition_timestamp,
	get_physiology_timestamp,
	get_number_of_samples,
	get_available_channels,
	get_active_channels,
	get_channel_mask,
	get_discard_pre,
	get_discard_post,
	get_center_sample,
	get_encoding_space_ref,
	get_trajectory_dimensions,
	get_sample_time_us,
	get_position,
	get_read_dir,
	get_phase_dir,
	get_slice_dir,
	get_patient_table_position,
	get_encoding_counters,
	get_user_ints,
	get_user_floats,
	get_acquisitionheader,
	is_flag_set,
	length,
	size,
	start,
	done,
	next,
	eltype,
	get_field_strength,
	get_num_rx_channels,
	get_enc_fov,
	get_enc_matrix,
	get_recon_fov,
	get_recon_matrix,
	get_encoding_limits,
	get_tr,
	get_te,
	get_ti,
	get_flip_angle,
	read_all_data;


# I should find a way to do this better
const libismrmrdwrap="/home/knowles/src/my_src/utils/io/ISMRMRD.jl/src/libismrmrdwrapper.so";

# static array lengths
const	ISMRMRD_USER_INTS = 8;
const	ISMRMRD_USER_FLOATS = 8;
const	ISMRMRD_PHYS_STAMPS = 3;
const	ISMRMRD_CHANNEL_MASKS = 16;
const	ISMRMRD_NDARRAY_MAXDIM = 7;
const	ISMRMRD_POSITION_LENGTH = 3;
const	ISMRMRD_DIRECTION_LENGTH = 3;

#flags ::Type{UInt64}
const   ISMRMRD_ACQ_FIRST_IN_ENCODE_STEP1               = 0x0000000000000001;
const   ISMRMRD_ACQ_LAST_IN_ENCODE_STEP1                = 0x0000000000000002;
const   ISMRMRD_ACQ_FIRST_IN_ENCODE_STEP2               = 0x0000000000000003;
const   ISMRMRD_ACQ_LAST_IN_ENCODE_STEP2                = 0x0000000000000004;
const   ISMRMRD_ACQ_FIRST_IN_AVERAGE                    = 0x0000000000000005;
const   ISMRMRD_ACQ_LAST_IN_AVERAGE                     = 0x0000000000000006;
const   ISMRMRD_ACQ_FIRST_IN_SLICE                      = 0x0000000000000007;
const   ISMRMRD_ACQ_LAST_IN_SLICE                       = 0x0000000000000008;
const   ISMRMRD_ACQ_FIRST_IN_CONTRAST                   = 0x0000000000000009;
const   ISMRMRD_ACQ_LAST_IN_CONTRAST                    = 0x000000000000000A;
const   ISMRMRD_ACQ_FIRST_IN_PHASE                      = 0x000000000000000B;
const   ISMRMRD_ACQ_LAST_IN_PHASE                       = 0x000000000000000C;
const   ISMRMRD_ACQ_FIRST_IN_REPETITION                 = 0x000000000000000D;
const   ISMRMRD_ACQ_LAST_IN_REPETITION                  = 0x000000000000000E;
const   ISMRMRD_ACQ_FIRST_IN_SET                        = 0x000000000000000F;
const   ISMRMRD_ACQ_LAST_IN_SET                         = 0x0000000000000010;
const   ISMRMRD_ACQ_FIRST_IN_SEGMENT                    = 0x0000000000000011;
const   ISMRMRD_ACQ_LAST_IN_SEGMENT                     = 0x0000000000000012;
const   ISMRMRD_ACQ_IS_NOISE_MEASUREMENT                = 0x0000000000000013;
const   ISMRMRD_ACQ_IS_PARALLEL_CALIBRATION             = 0x0000000000000014;
const   ISMRMRD_ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING = 0x0000000000000015;
const   ISMRMRD_ACQ_IS_REVERSE                          = 0x0000000000000016;
const   ISMRMRD_ACQ_IS_NAVIGATION_DATA                  = 0x0000000000000017;
const   ISMRMRD_ACQ_IS_PHASECORR_DATA                   = 0x0000000000000018;
const   ISMRMRD_ACQ_LAST_IN_MEASUREMENT                 = 0x0000000000000019;
const   ISMRMRD_ACQ_IS_HPFEEDBACK_DATA                  = 0x000000000000001A;
const   ISMRMRD_ACQ_IS_DUMMYSCAN_DATA                   = 0x000000000000001B;
const   ISMRMRD_ACQ_IS_RTFEEDBACK_DATA                  = 0x000000000000001C;
const   ISMRMRD_ACQ_IS_SURFACECOILCORRECTIONSCAN_DATA   = 0x000000000000001D;

const   ISMRMRD_ACQ_COMPRESSION1                        = 0x0000000000000035;
const   ISMRMRD_ACQ_COMPRESSION2                        = 0x0000000000000036;
const   ISMRMRD_ACQ_COMPRESSION3                        = 0x0000000000000037;
const   ISMRMRD_ACQ_COMPRESSION4                        = 0x0000000000000038;
const   ISMRMRD_ACQ_USER1                               = 0x0000000000000039;
const   ISMRMRD_ACQ_USER2                               = 0x000000000000003A;
const   ISMRMRD_ACQ_USER3                               = 0x000000000000003B;
const   ISMRMRD_ACQ_USER4                               = 0x000000000000003C;
const   ISMRMRD_ACQ_USER5                               = 0x000000000000003D;
const   ISMRMRD_ACQ_USER6                               = 0x000000000000003E;
const   ISMRMRD_ACQ_USER7                               = 0x000000000000003F;
const   ISMRMRD_ACQ_USER8                               = 0x0000000000000040;


# Error codes
const ISMRMRD_BEGINERROR 	=-1;
const ISMRMRD_NOERROR		= 0;
const ISMRMRD_MEMORYERROR	= 1;
const ISMRMRD_FILEERROR		= 2;
const ISMRMRD_TYPEERROR		= 3;
const ISMRMRD_RUNTIMEERROR 	= 4;
const ISMRMRD_HDF5ERROR 	= 5;
const ISMRMRD_ENDERROR		= 6;

#  ISMRMRD data types
const ISMRMRD_USHORT   = 1;# /**< corresponds to uint16_t */
const ISMRMRD_SHORT    = 2;# /**< corresponds to int16_t */
const ISMRMRD_UINT     = 3;# /**< corresponds to uint32_t */
const ISMRMRD_INT      = 4;# /**< corresponds to int32_t */
const ISMRMRD_FLOAT    = 5;# /**< corresponds to float */
const ISMRMRD_DOUBLE   = 6;# /**< corresponds to double */
const ISMRMRD_CXFLOAT  = 7;# /**< corresponds to complex float */
const ISMRMRD_CXDOUBLE = 8;# /**< corresponds to complex double */


# Image Types
const ISMRMRD_IMTYPE_MAGNITUDE = 1;
const ISMRMRD_IMTYPE_PHASE     = 2;
const ISMRMRD_IMTYPE_REAL      = 3;
const ISMRMRD_IMTYPE_IMAG      = 4;
const ISMRMRD_IMTYPE_COMPLEX   = 5;

# Image Flags
const ISMRMRD_IMAGE_IS_NAVIGATION_DATA =  1;
const ISMRMRD_IMAGE_USER1              = 57;
const ISMRMRD_IMAGE_USER2              = 58;
const ISMRMRD_IMAGE_USER3              = 59;
const ISMRMRD_IMAGE_USER4              = 60;
const ISMRMRD_IMAGE_USER5              = 61;
const ISMRMRD_IMAGE_USER6              = 62;
const ISMRMRD_IMAGE_USER7              = 63;
const ISMRMRD_IMAGE_USER8              = 64;



type ISMRMRD_dataset
	ptr::Ptr{Void};

	function ISMRMRD_dataset(filename::String,g::String="dataset",a::String="r")
	    if a == "w"
		create_file = 1;
	    else
		create_file = 0;
	    end


	    if (create_file == 0) && !isfile(filename)
	      error("Cannot open ", filename, " for reading. Does file exist?");
	    end

	    ptr = ccall((:jl_ismrmrd_create_dataset, libismrmrdwrap),
			Ptr{Void},(Cstring, Cstring, Cint), filename,g, create_file);

	    smart_ptr = new(ptr);
	    finalizer(smart_ptr,
	       obj -> ccall((:jl_ismrmrd_delete_dataset,libismrmrdwrap),
		    Void, (Ptr{Void},), obj.ptr));

	    return smart_ptr;
	end
end

type ISMRMRD_acquisition
	ptr::Ptr{Void};

	function ISMRMRD_acquisition()
		ptr = ccall((:jl_ismrmrd_create_acquisition, libismrmrdwrap),
	      		Ptr{Void},() );
		smart_ptr = new(ptr);
		finalizer(smart_ptr, obj -> 
	    		ccall((:jl_ismrmrd_delete_acquisition,libismrmrdwrap),
			Void, (Ptr{Void},), obj.ptr));
		return smart_ptr;
	end

	function ISMRMRD_acquisition(acq::ISMRMRD_acquisition)
		ptr = ccall((:jl_ismrmrd_copy_acquisition, libismrmrdwrap),
			Ptr{Void},(Ptr{Void},), acq.ptr);
		smart_ptr = new(ptr);
		finalizer(smart_ptr, obj -> 
	    		ccall((:jl_ismrmrd_delete_acquisition,libismrmrdwrap),
			Void, (Ptr{Void},), obj.ptr));
		return smart_ptr;
	end

	function ISMRMRD_acquisition(samples::Int32, 
			      channels::Int32, traj_dims::Int32)
		ptr = ccall((:jl_ismrmrd_setup_acquisition, libismrmrdwrap),
			Ptr{Void},(Cint, Cint, Cint), 
				samples, channels, traj_dims);
		smart_ptr = new(ptr);
		finalizer(smart_ptr, obj -> 
	    		ccall((:jl_ismrmrd_delete_acquisition,libismrmrdwrap),
			Void, (Ptr{Void},), obj.ptr));
		return smart_ptr;
	end

	function ISMRMRD_acquisition(ptr::Ptr{Void})
		return new(ptr);
	end

end

"""
ISMRMRD_image type
"""
type ISMRMRD_image
  ptr::Ptr{Void}
  function ISMRMRD_image()
    ptr = ccall((:jl_ismrmrd_create_image,libismrmrdwrap),
		Ptr{Void},());
    smart_ptr = new(ptr);
    finalizer(smart_ptr,obj -> ccall((:jl_ismrmrmd_delete_image,libismrmrdwrap),
				      Void, (Ptr{Void},),obj.ptr));
    return smart_ptr;
				     
  end
  
  function ISMRMRD_image(ptr::Ptr{Void})
    return new(ptr);
  end
end

type ISMRMRD_EncodingCounters
	kspace_encode_step_1::Int16
	kspace_encode_step_2::Int16
	average::Int16
	slice::Int16
	contrast::Int16
	phase::Int16
	repetition::Int16
	set::Int16
	segment::Int16
	user::Array{Int16,1};
end

type ISMRMRD_acquisitionheader
	version::Int64;
	flags::UInt64;
	measurement_uid::Int64;
	scan_counter::Int64;
	acquisition_time_stamp::Int64;
	physiology_time_stamp::Array{Int64,1};
	number_of_samples::Int64;
	available_channels::Int64;
	active_channels::Int64;
	channel_mask::Array{Int64,1};
	discard_pre::Int64;
	discard_post::Int64;
	center_sample::Int64;
	encoding_space_ref::Int64;
	trajectory_dimensions::Int64;
	sample_time_us::Float32;
	position::Array{Float32,1};
	read_dir::Array{Float32,1};
	phase_dir::Array{Float32,1};
	slice_dir::Array{Float32,1};
	patient_table_position::Array{Float32,1};
	idx::ISMRMRD_EncodingCounters;
	user_int::Array{Int32,1};
	user_float::Array{Int64,1};
end

#Methods:Dataset
"""
read_xml_header(dset::ISMRMRD_dataset)
- reads the ismrmrmd header file and returns the XML data as a string
"""
function read_xml_header(dset::ISMRMRD_dataset)
	buffer = Array{UInt8,1}(10000);

	sz = ccall((:jl_ismrmrd_read_header, libismrmrdwrap),
		Cint,
		(Ptr{Void},Ptr{UInt8}),
		dset.ptr,buffer);

	return unsafe_string(pointer(buffer),sz);
end

"""
write_xml_header(dset::ISMRMRD_dataset)
- writes an XML header to a ismrmrmd file
"""
function write_xml_header(dset::ISMRMRD_dataset, xml_hdr::String)
	ccall((:jl_ismrmrd_write_header, libismrmrdwrap),
		Void,
		(Ptr{Void},Cstring),
		dset.ptr,xml_hdr);
end

"""
append_dataset(dset::ISMRMRD_dataset, acq::ISMRMRD_acquisition)
- pushes an acquisition to the dataset
"""
function append_dataset(dset::ISMRMRD_dataset, acq::ISMRMRD_acquisition)
	ccall((:jl_ismrmrd_append_acquisition, libismrmrdwrap),
		Void,
		(Ptr{Void},Ptr{Void}),
		dset.ptr,acq.ptr);
end

"""
read_acquisition(dset::ISMRMRD_dataset, idx::Int32)
- returns an ISMRMRD_acquisition type from the dataset at index idx
- note that idx should be in the range of 1:nAcqs
"""
function read_acquisition(dset::ISMRMRD_dataset, idx::Int)
	#checks
	num_acqs = get_number_of_acquisitions(dset);
	if idx < 1 || idx > num_acqs
		error("idx $idx is not a valid index number");
	end

	#convert to c style as expected by libismrmrdwrap
	idx -= 1;

	ptr=ccall((:jl_ismrmrd_read_acquisition, libismrmrdwrap),
		Ptr{Void},
		(Ptr{Void},Cint),
		dset.ptr,idx);

	return ISMRMRD_acquisition(ptr);
end

"""
get_number_of_acquisitions(dset::ISMRMRD_dataset)
- returns the number of acquisitions in the dataset
"""
function get_number_of_acquisitions(dset::ISMRMRD_dataset)
	return ccall((:jl_ismrmrd_get_number_of_acquisitions, libismrmrdwrap),
		Cint,
		(Ptr{Void},),
		dset.ptr);
end

#
#methods: acquisition
#
function get_number_of_data_samples(acq::ISMRMRD_acquisition)
	return ccall((:jl_ismrmrd_get_number_of_data_elements, libismrmrdwrap),
		Cint,
		(Ptr{Void},),
		acq.ptr);
end

function get_number_of_trajectory_samples(acq::ISMRMRD_acquisition)
	return ccall((:jl_ismrmrd_get_number_of_traj_elements, libismrmrdwrap),
		Cint,
		(Ptr{Void},),
		acq.ptr);
end

function get_data(acq::ISMRMRD_acquisition)
	data_ptr=ccall((:jl_ismrmrd_get_data_pointer, libismrmrdwrap),
		Ptr{Complex{Float32}},
		(Ptr{Void},),
		acq.ptr);
	ns = get_number_of_data_samples(acq);

	data = [];
	try
		data = unsafe_wrap(Array{Complex{Float32},1},data_ptr,ns);
	catch
		error("failed to load data for $ns samples");
	end

	return data;
end

function get_trajectory(acq::ISMRMRD_acquisition)
	data_ptr=ccall((:jl_ismrmrd_get_traj_pointer, libismrmrdwrap),
		Ptr{Complex{Float32}},
		(Ptr{Void},),
		acq.ptr);
	nt = get_number_of_trajectory_samples(acq);

	traj = [];
	try
		traj = unsafe_wrap(Array{Float32,1},data_ptr,nt);
	catch
		error("failed to load for $nt samples");
	end

	return traj;
end

function get_version(acq::ISMRMRD_acquisition)
  x = ccall((:jl_ismrmrd_header_version, libismrmrdwrap), Culonglong,(Ptr{Void},),acq.ptr);
	return convert(Int64,x);
end

function get_flags(acq::ISMRMRD_acquisition)
	return ccall((:jl_ismrmrd_header_flags, libismrmrdwrap), Culonglong,(Ptr{Void},),acq.ptr);
end

function get_meas_uid(acq::ISMRMRD_acquisition)
	x = ccall((:jl_ismrmrd_header_measurement_uid, libismrmrdwrap),Cuint,(Ptr{Void},),acq.ptr);
	return convert(Int64,x);
end

function get_scan_counter(acq::ISMRMRD_acquisition)
	x = ccall((:jl_ismrmrd_header_scan_counter,libismrmrdwrap),Cuint,(Ptr{Void},),acq.ptr);
	return convert(Int64, x);
end

function get_acquisition_timestamp(acq::ISMRMRD_acquisition)
	x = ccall((:jl_ismrmrd_header_acquisition_time_stamp,libismrmrdwrap),Cuint,(Ptr{Void},),acq.ptr);
	return convert(Int64, x);
end

function get_physiology_timestamp(acq::ISMRMRD_acquisition)
	pts = Array{UInt32,1}(3);
	pts[1] = ccall((:jl_ismrmrd_header_physiology_time_stamp,libismrmrdwrap),Cuint,(Ptr{Void},Cint),acq.ptr,0);
	pts[2] = ccall((:jl_ismrmrd_header_physiology_time_stamp,libismrmrdwrap),Cuint,(Ptr{Void},Cint),acq.ptr,1);
	pts[3] = ccall((:jl_ismrmrd_header_physiology_time_stamp,libismrmrdwrap),Cuint,(Ptr{Void},Cint),acq.ptr,2);
	return convert(Array{Int64,1},pts);
end

function get_number_of_samples(acq::ISMRMRD_acquisition)
	x = ccall((:jl_ismrmrd_header_number_of_samples,libismrmrdwrap),Cushort,(Ptr{Void},),acq.ptr);
	return convert(Int64, x);
end

function get_available_channels(acq::ISMRMRD_acquisition)
	x = ccall((:jl_ismrmrd_header_available_channels,libismrmrdwrap),Cushort,(Ptr{Void},),acq.ptr);
	return convert(Int64, x);
end

function get_active_channels(acq::ISMRMRD_acquisition)
	x = ccall((:jl_ismrmrd_header_active_channels,libismrmrdwrap),Cushort,(Ptr{Void},),acq.ptr);
	return convert(Int64, x);
end

function get_channel_mask(acq::ISMRMRD_acquisition)
	cmask = Array{UInt64,1}(ISMRMRD_CHANNEL_MASKS);
	for c=1:ISMRMRD_CHANNEL_MASKS
		cmask[c] = ccall((:jl_ismrmrd_header_channel_mask,libismrmrdwrap),Cushort,(Ptr{Void},Cint),acq.ptr,c-1);
	end
	return convert(Array{Int64,1}, cmask);
end

function get_discard_pre(acq::ISMRMRD_acquisition)
	x = ccall((:jl_ismrmrd_header_discard_pre,libismrmrdwrap),Cushort,(Ptr{Void},),acq.ptr);
	return convert(Int64, x);
end

function get_discard_post(acq::ISMRMRD_acquisition)
	x = ccall((:jl_ismrmrd_header_discard_post,libismrmrdwrap),Cushort,(Ptr{Void},),acq.ptr);
	return convert(Int64, x);
end

function get_center_sample(acq::ISMRMRD_acquisition)
	x = ccall((:jl_ismrmrd_header_center_sample,libismrmrdwrap),Cushort,(Ptr{Void},),acq.ptr);
	return convert(Int64, x);
end

function get_encoding_space_ref(acq::ISMRMRD_acquisition)
	x = ccall((:jl_ismrmrd_header_encoding_space_ref,libismrmrdwrap),Cushort,(Ptr{Void},),acq.ptr);
	return convert(Int64, x);
end

function get_trajectory_dimensions(acq::ISMRMRD_acquisition)
	x = ccall((:jl_ismrmrd_header_trajectory_dimensions,libismrmrdwrap),Cushort,(Ptr{Void},),acq.ptr);
	return convert(Int64, x);
end

function get_sample_time_us(acq::ISMRMRD_acquisition)
	x = ccall((:jl_ismrmrd_header_sample_time_us,libismrmrdwrap),Cfloat,(Ptr{Void},),acq.ptr);
	return convert(Float64, x)
end

function get_position(acq::ISMRMRD_acquisition)
	x = Array{Float32,1}(ISMRMRD_POSITION_LENGTH);
	for c=1:ISMRMRD_POSITION_LENGTH
		x[c] = ccall((:jl_ismrmrd_header_position,libismrmrdwrap),Cfloat,(Ptr{Void},Cint),acq.ptr,c-1);
	end
	return x;
end

function get_read_dir(acq::ISMRMRD_acquisition)
	x = Array{Float32,1}(ISMRMRD_DIRECTION_LENGTH);
	for c=1:ISMRMRD_DIRECTION_LENGTH
		x[c] = ccall((:jl_ismrmrd_header_read_dir,libismrmrdwrap),Cfloat,(Ptr{Void},Cint),acq.ptr,c-1);
	end
	return x;
end

function get_phase_dir(acq::ISMRMRD_acquisition)
	x = Array{Float32,1}(ISMRMRD_DIRECTION_LENGTH);
	for c=1:ISMRMRD_DIRECTION_LENGTH
		x[c] = ccall((:jl_ismrmrd_header_phase_dir,libismrmrdwrap),Cfloat,(Ptr{Void},Cint),acq.ptr,c-1);
	end
	return x;
end

function get_slice_dir(acq::ISMRMRD_acquisition)
	x = Array{Float32,1}(ISMRMRD_DIRECTION_LENGTH);
	for c=1:ISMRMRD_DIRECTION_LENGTH
		x[c] = ccall((:jl_ismrmrd_header_slice_dir,libismrmrdwrap),Cfloat,(Ptr{Void},Cint),acq.ptr,c-1);
	end
	return x;
end

function get_patient_table_position(acq::ISMRMRD_acquisition)
	x = Array{Float32,1}(ISMRMRD_DIRECTION_LENGTH);
	for c=1:ISMRMRD_DIRECTION_LENGTH
		x[c] = ccall((:jl_ismrmrd_header_patient_table_position,libismrmrdwrap),Cfloat,(Ptr{Void},Cint),acq.ptr,c-1);
	end
	return x;
end

"""
  returns the encoding counter for a given acquisition type
  Note: encoding counter is converted to Julia style 1-index,
  from the C-style 0-idx
"""
function get_encoding_counters(acq::ISMRMRD_acquisition)
	k1 = ccall((:jl_ismrmrd_get_ec_kspace_encode_step_1,libismrmrdwrap),Cushort,(Ptr{Void},),acq.ptr);
	k2 = ccall((:jl_ismrmrd_get_ec_kspace_encode_step_2,libismrmrdwrap),Cushort,(Ptr{Void},),acq.ptr);
	a = ccall((:jl_ismrmrd_get_ec_average,libismrmrdwrap),Cushort,(Ptr{Void},),acq.ptr);
	sli = ccall((:jl_ismrmrd_get_ec_slice,libismrmrdwrap),Cushort,(Ptr{Void},),acq.ptr);
	co = ccall((:jl_ismrmrd_get_ec_contrast,libismrmrdwrap),Cushort,(Ptr{Void},),acq.ptr);
	p = ccall((:jl_ismrmrd_get_ec_phase,libismrmrdwrap),Cushort,(Ptr{Void},),acq.ptr);
	r = ccall((:jl_ismrmrd_get_ec_repetition,libismrmrdwrap),Cushort,(Ptr{Void},),acq.ptr);
	set = ccall((:jl_ismrmrd_get_ec_set,libismrmrdwrap),Cushort,(Ptr{Void},),acq.ptr);
	seg = ccall((:jl_ismrmrd_get_ec_segment,libismrmrdwrap),Cushort,(Ptr{Void},),acq.ptr);
	user = Array{Int16,1}(ISMRMRD_USER_INTS)
	for c=1:ISMRMRD_USER_INTS
		user[c] =  ccall((:jl_ismrmrd_get_ec_user,libismrmrdwrap),Cushort,(Ptr{Void},Cint),acq.ptr,c-1);
	end

	# convert to signed ints so we don't have to read hex
	return ISMRMRD_EncodingCounters(convert(Int16,k1+1),
					convert(Int16,k2)+1,
					convert(Int16,a+1),
					convert(Int16,sli+1),
					convert(Int16,co+1),
					convert(Int16,p+1),
					convert(Int16,r+1),
					convert(Int16,set+1),
					convert(Int16,seg+1),
					convert(Array{Int16,1},user.+1));
end

function get_user_ints(acq::ISMRMRD_acquisition)
	x = Array{Int32,1}(ISMRMRD_USER_INTS);
	for c=1:ISMRMRD_USER_INTS
		x[c] = ccall((:jl_ismrmrd_header_user_float,libismrmrdwrap),Cint,(Ptr{Void},Cint),acq.ptr,c-1);
	end
	return x;
end
function get_user_floats(acq::ISMRMRD_acquisition)
	x = Array{Float32,1}(ISMRMRD_USER_FLOATS);
	for c=1:ISMRMRD_USER_FLOATS
		x[c] = ccall((:jl_ismrmrd_header_user_float,libismrmrdwrap),Cfloat,(Ptr{Void},Cint),acq.ptr,c-1);
	end
	return x;
end

function get_acquisitionheader(acq::ISMRMRD_acquisition)
	version = get_version(acq);
	flags = get_flags(acq);
	measurement_uid = get_meas_uid(acq);
	scan_counter = get_scan_counter(acq);
	acquisition_time_stamp = get_acquisition_timestamp(acq);
	physiology_time_stamp = get_physiology_timestamp(acq);
	number_of_samples = get_number_of_samples(acq);
	available_channels = get_available_channels(acq);
	active_channels = get_active_channels(acq);
	channel_mask = get_channel_mask(acq);
	discard_pre = get_discard_pre(acq);
	discard_post = get_discard_post(acq);
	center_sample=get_center_sample(acq);
	encoding_space_ref = get_encoding_space_ref(acq);
	trajectory_dimensions = get_trajectory_dimensions(acq);
	sample_time_us=get_sample_time_us(acq);
	pos= get_position(acq);
	read_dir=get_read_dir(acq);
	phase_dir=get_phase_dir(acq);
	slice_dir=get_slice_dir(acq);
	patient_table_position= get_patient_table_position(acq);
	idx= get_encoding_counters(acq);
	user_int=get_user_ints(acq);
	user_float=get_user_floats(acq)

	return ISMRMRD_acquisitionheader(version,flags,measurement_uid,scan_counter,
	acquisition_time_stamp,physiology_time_stamp,number_of_samples,available_channels,
	active_channels, channel_mask, discard_pre, discard_post, center_sample,
	encoding_space_ref, trajectory_dimensions, sample_time_us, pos, read_dir,
	phase_dir, slice_dir, patient_table_position, idx, user_int, user_float);
end

#
# flags
#
function is_flag_set(acq::ISMRMRD_acquisition, flag_val::UInt64)
	flags = get_flags(acq);
	bitmask = UInt64(1) << (flag_val - 1);
	return (flags & bitmask) > 0;
end

function is_flag_set(acq::ISMRMRD_acquisition, flag_vals::Array{UInt64})
	b = broadcast((x)->is_flag_set(acq,x),flag_vals);
	return isempty(find((x)->x==true,b)) ? false : true;
end

#
# dataset iterator interface
#
function length(dset::ISMRMRD_dataset)
	return get_number_of_acquisitions(dset);
end

function size(dset::ISMRMRD_dataset)
	return (length(dset),)
end

function start(dset::ISMRMRD_dataset)
	return 1;
end

function done(dset::ISMRMRD_dataset, curr_idx::Int64)
	num_acqs = length(dset);
	return curr_idx > num_acqs;
end

function next(dset::ISMRMRD_dataset, curr_idx::Int64)
	acq=read_acquisition(dset,curr_idx)
	next_idx = curr_idx + 1;
	return (acq, next_idx);
end

function eltype(dset::ISMRMRD_dataset)
	return ISMRMRD_acquisition;
end

#
# regex searches of the XML header
# Warning: these may not be consistent with your header
#s
function get_field_strength(xmlhdr::String)
	pat = r"\<systemFieldStrength_T\>(\d\.?\d+)\<\/systemFieldStrength_T\>"s;

	m = match(pat,xmlhdr);

	#better way to check?
	if isa(m,Void)
		return [];
	else
		return float.(m.captures);
	end
end

function get_num_rx_channels(xmlhdr::String)
	pat = r"\<receiverChannels\>(\d+)\<\/receiverChannels\>"s;

	m = match(pat,xmlhdr);

	#better way to check?
	if isa(m,Void)
		return [];
	else
		return broadcast((x)->parse(Int64,x,10),m.captures);
	end
end

function get_enc_fov(xmlhdr::String)
	pat = r"\<encodedSpace\>.*\<fieldOfView_mm\>\n\t+\<x\>(\w+\.?\w+?)\<\/x\>\n\t+\<y\>(\w+\.?\w+)\<\/y\>\n\t+\<z\>(\w+\.?\w+)\<\/z\>\n\t+\<\/fieldOfView_mm\>.*\<\/encodedSpace\>"s;

	m = match(pat,xmlhdr);

	#better way to check?
	if isa(m,Void)
		return [];
	else
		return float.(m.captures);
	end
end


function get_enc_matrix(xmlhdr::String)
	pat = r"\<encodedSpace\>\n\t+\<matrixSize\>\n\t+\<x\>(\d+)\<\/x\>\n\t+\<y\>(\d+)\<\/y\>\n\t+\<z\>(\d+)\<\/z\>\n\t+\<\/matrixSize\>.*\<\/encodedSpace\>"s;

	m = match(pat,xmlhdr);

	#better way to check?
	if isa(m,Void)
		return [];
	else
		return broadcast((x)->parse(Int64,x,10),m.captures);
	end
end


function get_recon_fov(xmlhdr::String)
	pat = r"\<reconSpace\>.*\<fieldOfView_mm\>\n\t+\<x\>(\w+\.?\w+?)\<\/x\>\n\t+\<y\>(\w+\.?\w+)\<\/y\>\n\t+\<z\>(\w+\.?\w+)\<\/z\>\n\t+\<\/fieldOfView_mm\>.*\<\/reconSpace\>"s;

	m = match(pat,xmlhdr);

	#better way to check?
	if isa(m,Void)
		return [];
	else
		return float.(m.captures);
	end
end

function get_recon_matrix(xmlhdr::String)
	pat = r"\<reconSpace\>\n\t+\<matrixSize\>\n\t+\<x\>(\d+)\<\/x\>\n\t+\<y\>(\d+)\<\/y\>\n\t+\<z\>(\d+)\<\/z\>\n\t+\<\/matrixSize\>.*\<\/reconSpace\>"s;

	m = match(pat,xmlhdr);

	#better way to check?
	if isa(m,Void)
		return [];
	else
		return broadcast((x)->parse(Int64,x,10),m.captures);
	end
end

function get_encoding_limits(xmlhdr::String)
	# i hope you like horrible looking regex!
	pat=r"\<encodingLimits\>\n\t+\<kspace_encoding_step_1\>\n\t+\<minimum\>(\d+)\<\/minimum\>\n\t+\<maximum\>(\d+)\<\/maximum\>\n\t+\<center\>(\d+)\<\/center\>\n\t+\<\/kspace_encoding_step_1\>\n\t+\<kspace_encoding_step_2\>\n\t+\<minimum\>(\d+)\<\/minimum\>\n\t+\<maximum\>(\d+)\<\/maximum\>\n\t+\<center\>(\d+)\<\/center\>\n\t+\<\/kspace_encoding_step_2\>\n\t+\<average\>\n\t+\<minimum\>(\d+)\<\/minimum\>\n\t+\<maximum\>(\d+)\<\/maximum\>\n\t+\<center\>(\d+)\<\/center\>\n\t+\<\/average\>\n\t+\<slice\>\n\t+\<minimum\>(\d+)\<\/minimum\>\n\t+\<maximum\>(\d+)\<\/maximum\>\n\t+\<center\>(\d+)\<\/center\>\n\t+\<\/slice\>\n\t+\<contrast\>\n\t+\<minimum\>(\d+)\<\/minimum\>\n\t+\<maximum\>(\d+)\<\/maximum\>\n\t+\<center\>(\d+)\<\/center\>\n\t+\<\/contrast\>\n\t+\<phase\>\n\t+\<minimum\>(\d+)\<\/minimum\>\n\t+\<maximum\>(\d+)\<\/maximum\>\n\t+\<center\>(\d+)\<\/center\>\n\t+\<\/phase\>\n\t+\<repetition\>\n\t+\<minimum\>(\d+)\<\/minimum\>\n\t+\<maximum\>(\d+)\<\/maximum\>\n\t+\<center\>(\d+)\<\/center\>\n\t+\<\/repetition\>\n\t+\<set\>\n\t+\<minimum\>(\d+)\<\/minimum\>\n\t+\<maximum\>(\d+)\<\/maximum\>\n\t+\<center\>(\d+)\<\/center\>\n\t+\<\/set\>\n\t+\<segment\>\n\t+\<minimum\>(\d+)\<\/minimum\>\n\t+\<maximum\>(\d+)\<\/maximum\>\n\t+\<center\>(\d+)\<\/center\>\n\t+\<\/segment\>\n\t+\<\/encodingLimits\>"s;

	m = match(pat,xmlhdr);

	#better way to check?
	if isa(m,Void)
		return [];
	end

	# sort
	lims = broadcast((x)->parse(Int64,x,10),m.captures); #float.(m.captures);
	minLims = lims[1:3:end];
	maxLims = lims[2:3:end];
	cenLims = lims[3:3:end];



	return (minLims,maxLims,cenLims);

end

function get_tr(xmlhdr::String)
	pat = r"\<TR\>(\d\.?\d+)\<\/TR\>"s;

	m = match(pat,xmlhdr);

	#better way to check?
	if isa(m,Void)
		return [];
	else
		return float.(m.captures);
	end
end

function get_te(xmlhdr::String)
	pat = r"\<TE\>(\d\.?\d+)\<\/TE\>"s;

	m = match(pat,xmlhdr);

	#better way to check?
	if isa(m,Void)
		return [];
	else
		return float.(m.captures);
	end
end

function get_ti(xmlhdr::String)
	pat = r"\<TI\>(\d\.?\d+)\<\/TI\>"s;

	m = match(pat,xmlhdr);

	#better way to check?
	if isa(m,Void)
		return [];
	else
		return float.(m.captures);
	end
end

function get_flip_angle(xmlhdr::String)
	pat = r"\<flipAngle_deg\>(\d\.?\d+)\<\/flipAngle_deg\>"s;

	m = match(pat,xmlhdr);

	#better way to check?
	if isa(m,Void)
		return [];
	else
		return float.(m.captures);
	end
end

"""
 quick solution to read all data into an array
 returns a 11D array [nkx,nc,nky,nkz,nave,nsli,ncont,nph,nrep,nset,nseg]
"""
function read_all_data(dset::ISMRMRD_dataset, skip_flags::Array{UInt64,1})

	num_acqs = get_number_of_acquisitions(dset);
	#create an array based on encoding limits / matrix size
	hdr = read_xml_header(dset);
	enc_mat = get_enc_matrix(hdr);
	ns = enc_mat[1];
	nc = get_num_rx_channels(hdr);
	nc = nc[1];
	(lmin,lmax,lcen)=get_encoding_limits(hdr);

	#create array
	sz  = (	ns*nc,
		lmax[1]-lmin[1]+1,
		lmax[2]-lmin[2]+1,
		lmax[3]-lmin[3]+1,
		lmax[4]-lmin[4]+1,
		lmax[5]-lmin[5]+1,
		lmax[6]-lmin[6]+1,
		lmax[7]-lmin[7]+1,
		lmax[8]-lmin[8]+1,
		lmax[9]-lmin[9]+1);

	data = zeros(Complex{Float32},ns*nc,num_acqs);
	for acq in dset
		acqdata = get_data(acq);

		#skip certain flags
		if  !is_flag_set(acq,skip_flags)
			#idx = get_encoding_counters(acq);
			idx = get_scan_counter(acq);
			data[:,idx]=acqdata;
			#data[	:,
			#	idx.kspace_encode_step_1+1,
			#	idx.kspace_encode_step_2+1,
			#	idx.average+1,
			#	idx.slice+1,
			#	idx.contrast+1,
			#	idx.phase+1,
			#	idx.repetition+1,
			#	idx.set+1,
			#	idx.segment+1] = acqdata;
		end
	end

	#reshape & return
	# TODO: permute data? It is quite memory intensive...
	data=reshape(data,(ns,nc,sz[2],sz[3],sz[4],sz[5],sz[6],sz[7],sz[8],sz[9],sz[10]));
	return data;
end

"""
"""
function read_all_data(dset::ISMRMRD_dataset)
	# default flags to skip
	skip_flags = 	[ISMRMRD_ACQ_IS_NOISE_MEASUREMENT,
			ISMRMRD_ACQ_IS_PARALLEL_CALIBRATION,
			ISMRMRD_ACQ_IS_NAVIGATION_DATA,
			ISMRMRD_ACQ_IS_PHASECORR_DATA,
			ISMRMRD_ACQ_IS_DUMMYSCAN_DATA,
			ISMRMRD_ACQ_IS_RTFEEDBACK_DATA,
			ISMRMRD_ACQ_IS_SURFACECOILCORRECTIONSCAN_DATA];

	return read_all_data(dset,skip_flags);
end

"""
reads in ISMRMRD images from a dataset
"""
function read_image(dset::ISMRMRD_dataset,group::String, idx::Int)
	#checks
	# TODO, add checks
	#num_acqs = get_number_of_acquisitions(dset);
	#if idx < 1 || idx > num_acqs
	#	error("idx $idx is not a valid index number");
	#end

	ptr=ccall((:jl_ismrmrd_read_image, libismrmrdwrap),
		Ptr{Void},
		(Ptr{Void},Cstring,Cint),
		dset.ptr,group,idx-1);

	return ISMRMRD_image(ptr);
end


"""
returns the type of image data

for reference:
ISMRMRD_USHORT   = 1;# /**< corresponds to uint16_t */
ISMRMRD_SHORT    = 2;# /**< corresponds to int16_t */
ISMRMRD_UINT     = 3;# /**< corresponds to uint32_t */
ISMRMRD_INT      = 4;# /**< corresponds to int32_t */
ISMRMRD_FLOAT    = 5;# /**< corresponds to float */
ISMRMRD_DOUBLE   = 6;# /**< corresponds to double */
ISMRMRD_CXFLOAT  = 7;# /**< corresponds to complex float */
ISMRMRD_CXDOUBLE = 8;# /**< corresponds to complex double */
"""
function get_image_data_type(img::ISMRMRD_image)
  data_type = ccall((:jl_ismrmrd_get_image_data_type,libismrmrdwrap),
		    Cint,
		    (Ptr{Void},),
		    img.ptr);
  return data_type;
end

"""
returns image dimensions
"""
function get_image_dimensions(img::ISMRMRD_image)
  dims = ccall((:jl_ismrmrd_get_image_dims,libismrmrdwrap),
	      Ptr{Cushort},
	      (Ptr{Void},),
	     img.ptr);

  return convert.(Int,unsafe_wrap(Array,dims,3,false));
end


"""
returns the image data pointeed to by the ISMRMRD_image type
"""
function get_data(img::ISMRMRD_image)
  data_type = get_image_data_type(img);
  dims = get_image_dimensions(img);

  if data_type == 1
    ptr = ccall((:jl_ismrmrd_get_ushort_data,libismrmrdwrap),
		Ptr{Cushort},
		(Ptr{Void},),
		img.ptr);
  elseif data_type == 2
    ptr = ccall((:jl_ismrmrd_get_short_data,libismrmrdwrap),
		Ptr{Cshort},
		(Ptr{Void},),
		img.ptr);
  elseif data_type == 3
    ptr = ccall((:jl_ismrmrd_get_uint_data,libismrmrdwrap),
		Ptr{Cuint},
		(Ptr{Void},),
		img.ptr);
  elseif data_type == 4
    ptr = ccall((:jl_ismrmrd_get_int_data,libismrmrdwrap),
		Ptr{Cint},
		(Ptr{Void},),
		img.ptr);
  elseif data_type == 5
    ptr = ccall((:jl_ismrmrd_get_float_data,libismrmrdwrap),
		Ptr{Cfloat},
		(Ptr{Void},),
		img.ptr);
  elseif data_type == 6
    ptr = ccall((:jl_ismrmrd_get_double_data,libismrmrdwrap),
		Ptr{Cdouble},
		(Ptr{Void},),
		img.ptr);
  elseif data_type == 7
    ptr = ccall((:jl_ismrmrd_get_complex_float_data,libismrmrdwrap),
		Ptr{Complex64},
		(Ptr{Void},),
		img.ptr);
  elseif data_type == 8
    ptr = ccall((:jl_ismrmrd_get_complex_double_data,libismrmrdwrap),
		Ptr{Complex128},
		(Ptr{Void},),
		img.ptr);
  else
    error("unknown data type");
  end

  # will return an array of type T, given a pointer of Ptr{T}
  return unsafe_wrap(Array,ptr,(dims[1],dims[2],dims[3]),true); # let julia take ownership of the memory

end


end #module
