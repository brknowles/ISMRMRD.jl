module ISMRMRD2

const libismrmrd="/home/knowles/src/my_src/utils/io/ISMRMRD2.jl/src/libismrmrdwrapper.so";

function __init__()
	if !isfile(libismrmrd)
		error("Could not find $libismrmrd");
	end
end

const ISMRMRD_USER_INTS = 8;
const ISMRMRD_USER_FLOATS = 8;
const ISMRMRD_PHYS_STAMPS = 3;
const ISMRMRD_CHANNEL_MASKS = 16;
const ISMRMRD_NDARRAY_MAXDIM = 7;
const ISMRMRD_POSITION_LENGTH = 3;
const ISMRMRD_DIRECTION_LENGTH = 3;

type ISMRMRD_dataset
	ptr::Ptr{Void};

	function ISMRMRD_dataset(filename::String,a::String="r")
	    if a == "w"
		create_file = 1;
	    else
		create_file = 0;
	    end

	    #TODO, check for file existence

	    ptr = ccall((:jl_ismrmrd_create_dataset, libismrmrd),
			Ptr{Void},(Cstring, Cint), filename, create_file);

	    smart_ptr = new(ptr);
	    finalizer(smart_ptr,
	       obj -> ccall((:jl_ismrmrd_delete_dataset,libismrmrd),
		    Void, (Ptr{Void},), obj.ptr));

	    return smart_ptr;
	end
end

type ISMRMRD_acquisition
	ptr::Ptr{Void};

	function ISMRMRD_acquisition()
		ptr = ccall((:jl_ismrmrd_create_acquisition, libismrmrd),
	      		Ptr{Void},() );
		smart_ptr = new(ptr);
		finalizer(smart_ptr, obj -> ccall((:jl_ismrmrd_delete_acquisition,libismrmrd),
			Void, (Ptr{Void},), obj.ptr));
		return smart_ptr;
	end

	function ISMRMRD_acquisition(acq::ISMRMRD_acquisition)
		ptr = ccall((:jl_ismrmrd_copy_acquisition, libismrmrd),
			Ptr{Void},(Ptr{Void},), acq.ptr);
		smart_ptr = new(ptr);
		finalizer(smart_ptr, obj -> ccall((:jl_ismrmrd_delete_acquisition,libismrmrd),
			Void, (Ptr{Void},), obj.ptr));
		return smart_ptr;
	end

	function ISMRMRD_acquisition(samples::Int32, channels::Int32, traj_dims::Int32)
		ptr = ccall((:jl_ismrmrd_setup_acquisition, libismrmrd),
			Ptr{Void},(Cint, Cint, Cint), samples, channels, traj_dims);
		smart_ptr = new(ptr);
		finalizer(smart_ptr, obj -> ccall((:jl_ismrmrd_delete_acquisition,libismrmrd),
			Void, (Ptr{Void},), obj.ptr));
		return smart_ptr;
	end

	function ISMRMRD_acquisition(ptr::Ptr{Void})
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

	sz = ccall((:jl_ismrmrd_read_header, libismrmrd),
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
	ccall((:jl_ismrmrd_write_header, libismrmrd),
		Void,
		(Ptr{Void},Cstring),
		dset.ptr,xml_hdr);
end

"""
append_dataset(dset::ISMRMRD_dataset, acq::ISMRMRD_acquisition)
- pushes an acquisition to the dataset
"""
function append_dataset(dset::ISMRMRD_dataset, acq::ISMRMRD_acquisition)
	ccall((:jl_ismrmrd_append_acquisition, libismrmrd),
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

	#convert to c style as expected by libismrmrd
	idx -= 1;

	ptr=ccall((:jl_ismrmrd_read_acquisition, libismrmrd),
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
	return ccall((:jl_ismrmrd_get_number_of_acquisitions, libismrmrd),
		Cint,
		(Ptr{Void},),
		dset.ptr);
end

#
#methods: acquisition
#
function get_number_of_data_samples(acq::ISMRMRD_acquisition)
	return ccall((:jl_ismrmrd_get_number_of_data_elements, libismrmrd),
		Cint,
		(Ptr{Void},),
		acq.ptr);
end

function get_number_of_trajectory_samples(acq::ISMRMRD_acquisition)
	return ccall((:jl_ismrmrd_get_number_of_traj_elements, libismrmrd),
		Cint,
		(Ptr{Void},),
		acq.ptr);
end

function get_data(acq::ISMRMRD_acquisition)
	data_ptr=ccall((:jl_ismrmrd_get_data_pointer, libismrmrd),
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

	#copy data out of the realm of C
	return copy(data);
end

function get_trajectory(acq::ISMRMRD_acquisition)
	data_ptr=ccall((:jl_ismrmrd_get_traj_pointer, libismrmrd),
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

	#copy data out of the realm of C
	return copy(traj);
end

function get_version(acq::ISMRMRD_acquisition)
  x = ccall((:jl_ismrmrd_header_version, libismrmrd), Culonglong,(Ptr{Void},),acq.ptr);
	return convert(Int64,x);
end

function get_flags(acq::ISMRMRD_acquisition)
	return ccall((:jl_ismrmrd_header_flags, libismrmrd), Culonglong,(Ptr{Void},),acq.ptr);
end

function get_meas_uid(acq::ISMRMRD_acquisition)
	x = ccall((:jl_ismrmrd_header_measurement_uid, libismrmrd),Cuint,(Ptr{Void},),acq.ptr);
	return convert(Int64,x);
end

function get_scan_counter(acq::ISMRMRD_acquisition)
	x = ccall((:jl_ismrmrd_header_scan_counter,libismrmrd),Cuint,(Ptr{Void},),acq.ptr);
	return convert(Int64, x);
end

function get_acquisition_timestamp(acq::ISMRMRD_acquisition)
	x = ccall((:jl_ismrmrd_header_acquisition_time_stamp,libismrmrd),Cuint,(Ptr{Void},),acq.ptr);
	return convert(Int64, x);
end

function get_physiology_timestamp(acq::ISMRMRD_acquisition)
	pts = Array{UInt32,1}(3);
	pts[1] = ccall((:jl_ismrmrd_header_physiology_time_stamp,libismrmrd),Cuint,(Ptr{Void},Cint),acq.ptr,0);
	pts[2] = ccall((:jl_ismrmrd_header_physiology_time_stamp,libismrmrd),Cuint,(Ptr{Void},Cint),acq.ptr,1);
	pts[3] = ccall((:jl_ismrmrd_header_physiology_time_stamp,libismrmrd),Cuint,(Ptr{Void},Cint),acq.ptr,2);
	return convert(Array{Int64,1},pts);
end

function get_number_of_samples(acq::ISMRMRD_acquisition)
	x = ccall((:jl_ismrmrd_header_number_of_samples,libismrmrd),Cushort,(Ptr{Void},),acq.ptr);
	return convert(Int64, x);
end

function get_available_channels(acq::ISMRMRD_acquisition)
	x = ccall((:jl_ismrmrd_header_available_channels,libismrmrd),Cushort,(Ptr{Void},),acq.ptr);
	return convert(Int64, x);
end

function get_active_channels(acq::ISMRMRD_acquisition)
	x = ccall((:jl_ismrmrd_header_active_channels,libismrmrd),Cushort,(Ptr{Void},),acq.ptr);
	return convert(Int64, x);
end

function get_channel_mask(acq::ISMRMRD_acquisition)
	cmask = Array{UInt64,1}(ISMRMRD_CHANNEL_MASKS);
	for c=1:ISMRMRD_CHANNEL_MASKS
		cmask[c] = ccall((:jl_ismrmrd_header_channel_mask,libismrmrd),Cushort,(Ptr{Void},Cint),acq.ptr,c-1);
	end
	return convert(Array{Int64,1}, cmask);
end

function get_discard_pre(acq::ISMRMRD_acquisition)
	x = ccall((:jl_ismrmrd_header_discard_pre,libismrmrd),Cushort,(Ptr{Void},),acq.ptr);
	return convert(Int64, x);
end

function get_discard_post(acq::ISMRMRD_acquisition)
	x = ccall((:jl_ismrmrd_header_discard_post,libismrmrd),Cushort,(Ptr{Void},),acq.ptr);
	return convert(Int64, x);
end

function get_center_sample(acq::ISMRMRD_acquisition)
	x = ccall((:jl_ismrmrd_header_center_sample,libismrmrd),Cushort,(Ptr{Void},),acq.ptr);
	return convert(Int64, x);
end

function get_encoding_space_ref(acq::ISMRMRD_acquisition)
	x = ccall((:jl_ismrmrd_header_encoding_space_ref,libismrmrd),Cushort,(Ptr{Void},),acq.ptr);
	return convert(Int64, x);
end

function get_trajectory_dimensions(acq::ISMRMRD_acquisition)
	x = ccall((:jl_ismrmrd_header_trajectory_dimensions,libismrmrd),Cushort,(Ptr{Void},),acq.ptr);
	return convert(Int64, x);
end

function get_sample_time_us(acq::ISMRMRD_acquisition)
	x = ccall((:jl_ismrmrd_header_sample_time_us,libismrmrd),Cfloat,(Ptr{Void},),acq.ptr);
	return convert(Float64, x)
end

function get_position(acq::ISMRMRD_acquisition)
	x = Array{Float32,1}(ISMRMRD_POSITION_LENGTH);
	for c=1:ISMRMRD_POSITION_LENGTH
		x[c] = ccall((:jl_ismrmrd_header_position,libismrmrd),Cfloat,(Ptr{Void},Cint),acq.ptr,c-1);
	end
	return x;
end

function get_read_dir(acq::ISMRMRD_acquisition)
	x = Array{Float32,1}(ISMRMRD_DIRECTION_LENGTH);
	for c=1:ISMRMRD_DIRECTION_LENGTH
		x[c] = ccall((:jl_ismrmrd_header_read_dir,libismrmrd),Cfloat,(Ptr{Void},Cint),acq.ptr,c-1);
	end
	return x;
end

function get_phase_dir(acq::ISMRMRD_acquisition)
	x = Array{Float32,1}(ISMRMRD_DIRECTION_LENGTH);
	for c=1:ISMRMRD_DIRECTION_LENGTH
		x[c] = ccall((:jl_ismrmrd_header_phase_dir,libismrmrd),Cfloat,(Ptr{Void},Cint),acq.ptr,c-1);
	end
	return x;
end

function get_slice_dir(acq::ISMRMRD_acquisition)
	x = Array{Float32,1}(ISMRMRD_DIRECTION_LENGTH);
	for c=1:ISMRMRD_DIRECTION_LENGTH
		x[c] = ccall((:jl_ismrmrd_header_slice_dir,libismrmrd),Cfloat,(Ptr{Void},Cint),acq.ptr,c-1);
	end
	return x;
end

function get_patient_table_position(acq::ISMRMRD_acquisition)
	x = Array{Float32,1}(ISMRMRD_DIRECTION_LENGTH);
	for c=1:ISMRMRD_DIRECTION_LENGTH
		x[c] = ccall((:jl_ismrmrd_header_patient_table_position,libismrmrd),Cfloat,(Ptr{Void},Cint),acq.ptr,c-1);
	end
	return x;
end

function get_encoding_counters(acq::ISMRMRD_acquisition)
	k1 = ccall((:jl_ismrmrd_get_ec_kspace_encode_step_1,libismrmrd),Cushort,(Ptr{Void},),acq.ptr);
	k2 = ccall((:jl_ismrmrd_get_ec_kspace_encode_step_2,libismrmrd),Cushort,(Ptr{Void},),acq.ptr);
	a = ccall((:jl_ismrmrd_get_ec_average,libismrmrd),Cushort,(Ptr{Void},),acq.ptr);
	sli = ccall((:jl_ismrmrd_get_ec_slice,libismrmrd),Cushort,(Ptr{Void},),acq.ptr);
	c = ccall((:jl_ismrmrd_get_ec_contrast,libismrmrd),Cushort,(Ptr{Void},),acq.ptr);
	p = ccall((:jl_ismrmrd_get_ec_phase,libismrmrd),Cushort,(Ptr{Void},),acq.ptr);
	r = ccall((:jl_ismrmrd_get_ec_repetition,libismrmrd),Cushort,(Ptr{Void},),acq.ptr);
	set = ccall((:jl_ismrmrd_get_ec_set,libismrmrd),Cushort,(Ptr{Void},),acq.ptr);
	seg = ccall((:jl_ismrmrd_get_ec_segment,libismrmrd),Cushort,(Ptr{Void},),acq.ptr);
	user = Array{Int16,1}(ISMRMRD_USER_INTS)
	for c=1:ISMRMRD_USER_INTS
		user[c] =  ccall((:jl_ismrmrd_get_ec_user,libismrmrd),Cushort,(Ptr{Void},Cint),acq.ptr,c-1);
	end

	# convert to signed ints so we don't have to read hex
	return ISMRMRD_EncodingCounters(convert(Int16,k1),
																	convert(Int16,k2),
																	convert(Int16,a),
																	convert(Int16,sli),
																	convert(Int16,c),
																	convert(Int16,p),
																	convert(Int16,r),
																	convert(Int16,set),
																	convert(Int16,seg),
																	convert(Array{Int16,1},user));
end

function get_user_ints(acq::ISMRMRD_acquisition)
	x = Array{Int32,1}(ISMRMRD_USER_INTS);
	for c=1:ISMRMRD_USER_INTS
		x[c] = ccall((:jl_ismrmrd_header_user_float,libismrmrd),Cint,(Ptr{Void},Cint),acq.ptr,c-1);
	end
	return x;
end
function get_user_floats(acq::ISMRMRD_acquisition)
	x = Array{Float32,1}(ISMRMRD_USER_FLOATS);
	for c=1:ISMRMRD_USER_FLOATS
		x[c] = ccall((:jl_ismrmrd_header_user_float,libismrmrd),Cfloat,(Ptr{Void},Cint),acq.ptr,c-1);
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


end #module
