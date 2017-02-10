#ifndef _GDCM_WRAPPER_
#define _GDCM_WRAPPER_ 1

#include "ismrmrd/ismrmrd.h"

#ifdef __cplusplus
extern "C" {
#endif
	/*
	 *
	 *	ISMRMRD::Dataset interface
	 *
	*/
	void*  jl_ismrmrd_create_dataset(const char*, int );
	void   jl_ismrmrd_delete_dataset(void* );
	/*int jl_ismrmrd_is_dataset_open(void *);*/
  int jl_ismrmrd_read_header(void*, char* );
	/*void* jl_ismrmrd_deserialize_header(const char* );*/
	void jl_ismrmrd_write_header(void*, const char* );
	void jl_ismrmrd_append_acquisition(void*, void *);
	void* jl_ismrmrd_read_acquisition(void*, int32_t);
	int jl_ismrmrd_get_number_of_acquisitions(void* );

	/*
	 *
	 *	ISMRMRD::Acquisition interface
	 *
	 */
	void* jl_ismrmrd_create_acquisition(void );
	void* jl_ismrmrd_setup_acquisition(int ,int, int);
	void* jl_ismrmrd_copy_acquisition_copy(void* );
	void jl_ismrmrd_delete_acquisition(void *);

	/*
	 * Accessors
	 */
	uint64_t jl_ismrmrd_header_version(void*);
	uint64_t jl_ismrmrd_header_flags(void*);
	uint32_t jl_ismrmrd_header_measurement_uid(void*);
	uint32_t jl_ismrmrd_header_scan_counter(void*);
	uint32_t jl_ismrmrd_header_acquisition_time_stamp(void*);
	uint32_t jl_ismrmrd_header_physiology_time_stamp(void*, int );
	const uint16_t  jl_ismrmrd_header_number_of_samples(void*);
	uint16_t jl_ismrmrd_header_available_channels(void*);
	const uint16_t jl_ismrmrd_header_active_channels(void*);
	const uint64_t jl_ismrmrd_header_channel_mask(void*, int );
	uint16_t jl_ismrmrd_header_discard_pre(void*);
	uint16_t jl_ismrmrd_header_discard_post(void*);
	uint16_t jl_ismrmrd_header_center_sample(void*);
	uint16_t jl_ismrmrd_header_encoding_space_ref(void*);
	const uint16_t jl_ismrmrd_header_trajectory_dimensions(void*);
	float jl_ismrmrd_header_sample_time_us(void*);
	float jl_ismrmrd_header_position(void*, int );
	float jl_ismrmrd_header_read_dir(void*, int );
	float jl_ismrmrd_header_phase_dir(void*, int );
	float jl_ismrmrd_header_slice_dir(void*, int );
	float jl_ismrmrd_header_patient_table_position(void*, int );
	/*void* jl_ismrmrd_header_idx(void*);*/
	int jl_ismrmrd_header_user_int(void*, int );
	float jl_ismrmrd_header_user_float(void*, int );

	/*
	 * EncodingCounter
	 */
	 uint16_t jl_ismrmrd_get_ec_kspace_encode_step_1(void *);
	 uint16_t jl_ismrmrd_get_ec_kspace_encode_step_2(void *);
	 uint16_t jl_ismrmrd_get_ec_average(void *);
	 uint16_t jl_ismrmrd_get_ec_slice(void *);
	 uint16_t jl_ismrmrd_get_ec_contrast(void *);
	 uint16_t jl_ismrmrd_get_ec_phase(void *);
	 uint16_t jl_ismrmrd_get_ec_repetition(void *);
	 uint16_t jl_ismrmrd_get_ec_set(void *);
	 uint16_t jl_ismrmrd_get_ec_segment(void *);
	 uint16_t jl_ismrmrd_get_ec_user(void *, int);

	/*
	 * sizes
	 */
	void jl_ismrmrd_resize(void*, uint16_t, uint16_t, uint16_t);
	size_t jl_ismrmrd_get_number_of_data_elements(void* );
	size_t jl_ismrmrd_get_number_of_traj_elements(void* );
	size_t jl_ismrmrd_get_data_size(void *);
	size_t jl_ismrmrd_get_traj_size(void *);

	/*
	 * data/traj accessors
	 */
	/*
	ISMRMRD_AcquisitionHeader jl_ismrmrd_get_head(void *);
	void jl_ismrmrd_set_head(void*, ISMRMRD_AcquisitionHeader );
	*/

	complex_float_t* jl_ismrmrd_get_data_pointer(void* );
	float* jl_ismrmrd_get_traj_pointer(void* );

	/*
	 * flags
	 */
	int jl_ismrmrd_is_flag_set(void*, const uint16_t );
	void jl_ismrmrd_set_flag(void*, const uint16_t );
	void jl_ismrmrd_clear_flag(void*, const uint16_t);
	void jl_ismrmrd_clear_all_flags(void*);

	/*
	 * channels
	 */
	int jl_ismrmrd_is_channel_active(void*, uint16_t );
	void jl_ismrmrd_set_channel_active(void*, uint16_t );
	void jl_ismrmrd_set_channel_not_active(void*, uint16_t);
	void jl_ismrmrd_set_all_channels_not_active(void*);



#ifdef __cplusplus
}
#endif

#endif
