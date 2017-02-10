#include "ismrmrd/ismrmrd.h"
#include "ismrmrd/dataset.h"
#include "ismrmrd/xml.h"
#include <iostream>
#include <string>
#include "ismrmrd_wrapper.h"

//#define DEBUG 1

// Reader wrappers
void*
jl_ismrmrd_create_dataset(const char * filename,int create_file) {
	bool bCreate = (create_file > 0);

#ifdef DEBUG
	std::cout << "Creating ismrmrd_dataset \n";
#endif

	/*
	 * Dataset is initialied and file is read
	 */
	try{
		ISMRMRD::Dataset* dset =
			new ISMRMRD::Dataset(filename,"dataset",bCreate);
		return reinterpret_cast<void*>(dset);
	} catch (...) {
		std::cout<<" jl_ismrmrd_create_dataset: " <<
			"Error opening dataset file " <<
			filename <<std::endl;
		return NULL;
	}
}

void
jl_ismrmrd_delete_dataset(void* v_ds) {
	ISMRMRD::Dataset* dset = reinterpret_cast<ISMRMRD::Dataset*>(v_ds);

#ifdef DEBUG
	std::cout << "do we delete ismrmrd_dataset? \n";
#endif
	if (dset) {
#ifdef DEBUG
	std::cout << "deleting ismrmrd_dataset \n";
#endif
		delete dset;
	}
	return;
}

int
jl_ismrmrd_read_header(void* v_ds, char* hdr) {
	ISMRMRD::Dataset* dset = reinterpret_cast<ISMRMRD::Dataset*>(v_ds);
	std::string shdr;
	dset->readHeader(shdr);

	//copy...
	for(int i=0; i<(int) shdr.size(); i++) {
		hdr[i] = shdr[i];
	}

	return shdr.size();
}

void
jl_ismrmrd_write_header(void* v_ds, const char* xml_hdr) {
	ISMRMRD::Dataset* dset = reinterpret_cast<ISMRMRD::Dataset*>(v_ds);
	std::string shdr(xml_hdr);
	dset->writeHeader(shdr);

	return;
}
/*
void*
jl_ismrmrd_deserialize_header(const char* xml_header) {
	ISMRMRD::IsmrmrdHeader structHead;
	ISMRMRD::deserialize(xml_header,structHead);
	return reinterpret_cast<void*>(&structHead);
}
*/

void
jl_ismrmrd_append_acquisition(void* v_ds, void* v_aq) {
	if(!v_ds || !v_aq)
		return;

	ISMRMRD::Dataset* dset    = reinterpret_cast<ISMRMRD::Dataset*>(v_ds);
	ISMRMRD::Acquisition* acq = reinterpret_cast<ISMRMRD::Acquisition*>(v_aq);
	dset->appendAcquisition(*acq);

	return;

}

void* jl_ismrmrd_read_acquisition(void* v_ds, int32_t idx) {
	if(!v_ds)
		return NULL;

	ISMRMRD::Dataset* dset	= reinterpret_cast<ISMRMRD::Dataset*>(v_ds);
	ISMRMRD::Acquisition* acq = new ISMRMRD::Acquisition();
	dset->readAcquisition((uint32_t)idx, *acq);

	return acq;
}

int
jl_ismrmrd_get_number_of_acquisitions(void* v_ds ) {
	if(!v_ds)
		return 1;

	ISMRMRD::Dataset* dset = reinterpret_cast<ISMRMRD::Dataset*>(v_ds);

	return (int)dset->getNumberOfAcquisitions();
}

//acquisition
void*
jl_ismrmrd_create_acquisition(void) {
	ISMRMRD::Acquisition* aq = new ISMRMRD::Acquisition();
	return reinterpret_cast<void*>(aq);
}

void*
jl_ismrmrd_copy_acquisition(void * v_aq) {
	if(!v_aq)
		return NULL;

	ISMRMRD::Acquisition* acqsrc = reinterpret_cast<ISMRMRD::Acquisition*>(v_aq);
	ISMRMRD::Acquisition* acqdest = new ISMRMRD::Acquisition(*acqsrc);
	return reinterpret_cast<void*>(acqdest);
}

void*
jl_ismrmrd_setup_acquisition(int num_samples, int active_channels, int trajectory_dimensions) {

	ISMRMRD::Acquisition* acq = new ISMRMRD::Acquisition(	(uint16_t)num_samples,
								(uint16_t)active_channels,
								(uint16_t)trajectory_dimensions);
	return reinterpret_cast<void*>(acq);
}

void
jl_ismrmrd_delete_acquisition(void* v_aq) {
	if(!v_aq){
		return;
	}

	ISMRMRD::Acquisition* aq = reinterpret_cast<ISMRMRD::Acquisition*>(v_aq);
	delete aq;

	return;
}

//header
uint64_t
jl_ismrmrd_header_version(void* v_acq) {
	if(!v_acq){
		return 0;
	}
	ISMRMRD::Acquisition* acq = reinterpret_cast<ISMRMRD::Acquisition*>(v_acq);
	return acq->version();
}

uint64_t
jl_ismrmrd_header_flags(void* v_acq) {
	if(!v_acq){
		return 0;
	}
	ISMRMRD::Acquisition* acq = reinterpret_cast<ISMRMRD::Acquisition*>(v_acq);
	return acq->flags();
}

uint32_t
jl_ismrmrd_header_measurement_uid(void* v_acq) {
	if(!v_acq){
		return 0;
	}
	ISMRMRD::Acquisition* acq = reinterpret_cast<ISMRMRD::Acquisition*>(v_acq);
	return acq->measurement_uid();
}

uint32_t
jl_ismrmrd_header_scan_counter(void* v_acq) {
	if(!v_acq){
		return 0;
	}
	ISMRMRD::Acquisition* acq = reinterpret_cast<ISMRMRD::Acquisition*>(v_acq);
	return (int)acq->scan_counter();
}

uint32_t
jl_ismrmrd_header_acquisition_time_stamp(void* v_acq) {
	if(!v_acq){
		return 0;
	}
	ISMRMRD::Acquisition* acq = reinterpret_cast<ISMRMRD::Acquisition*>(v_acq);
	return (int)acq->acquisition_time_stamp();
}

uint32_t
jl_ismrmrd_header_physiology_time_stamp(void* v_acq, int idx) {
	if(!v_acq){
		return 0;
	}
	ISMRMRD::Acquisition* acq = reinterpret_cast<ISMRMRD::Acquisition*>(v_acq);
	return acq->physiology_time_stamp()[idx];
}

const uint16_t
jl_ismrmrd_header_number_of_samples(void* v_acq) {
	if(!v_acq){
		return 0;
	}
	ISMRMRD::Acquisition* acq = reinterpret_cast<ISMRMRD::Acquisition*>(v_acq);
	return acq->number_of_samples();
}

uint16_t
jl_ismrmrd_header_available_channels(void* v_acq) {
	if(!v_acq){
		return 0;
	}
	ISMRMRD::Acquisition* acq = reinterpret_cast<ISMRMRD::Acquisition*>(v_acq);
	return acq->available_channels();
}

const uint16_t
jl_ismrmrd_header_active_channels(void* v_acq) {
	if(!v_acq){
		return 0;
	}
	ISMRMRD::Acquisition* acq = reinterpret_cast<ISMRMRD::Acquisition*>(v_acq);
	return acq->active_channels();
}

const uint64_t
jl_ismrmrd_header_channel_mask(void* v_acq, int idx) {
	if(!v_acq){
		return 0;
	}
	ISMRMRD::Acquisition* acq = reinterpret_cast<ISMRMRD::Acquisition*>(v_acq);
	return acq->channel_mask()[idx];
}

uint16_t
jl_ismrmrd_header_discard_pre(void* v_acq) {
	if(!v_acq){
		return 0;
	}
	ISMRMRD::Acquisition* acq = reinterpret_cast<ISMRMRD::Acquisition*>(v_acq);
	return acq->discard_pre();
}

uint16_t
jl_ismrmrd_header_discard_post(void* v_acq) {
	if(!v_acq){
		return 0;
	}
	ISMRMRD::Acquisition* acq = reinterpret_cast<ISMRMRD::Acquisition*>(v_acq);
	return acq->discard_post();
}

uint16_t
jl_ismrmrd_header_center_sample(void* v_acq) {
	if(!v_acq){
		return 0;
	}
	ISMRMRD::Acquisition* acq = reinterpret_cast<ISMRMRD::Acquisition*>(v_acq);
	return acq->center_sample();
}

uint16_t
jl_ismrmrd_header_encoding_space_ref(void* v_acq) {
	if(!v_acq){
		return 0;
	}
	ISMRMRD::Acquisition* acq = reinterpret_cast<ISMRMRD::Acquisition*>(v_acq);
	return acq->encoding_space_ref();
}

const uint16_t
jl_ismrmrd_header_trajectory_dimensions(void* v_acq) {
	if(!v_acq){
		return 0;
	}
	ISMRMRD::Acquisition* acq = reinterpret_cast<ISMRMRD::Acquisition*>(v_acq);
	return acq->trajectory_dimensions();
}

float
jl_ismrmrd_header_sample_time_us(void* v_acq) {
	if(!v_acq){
		return 0.;
	}
	ISMRMRD::Acquisition* acq = reinterpret_cast<ISMRMRD::Acquisition*>(v_acq);
	return acq->sample_time_us();
}

float
jl_ismrmrd_header_position(void* v_acq, int idx) {
	if(!v_acq){
		return 0.;
	}
	ISMRMRD::Acquisition* acq = reinterpret_cast<ISMRMRD::Acquisition*>(v_acq);
	return acq->position()[idx];
}

float
jl_ismrmrd_header_read_dir(void* v_acq, int idx) {
	if(!v_acq){
		return 0.;
	}
	ISMRMRD::Acquisition* acq = reinterpret_cast<ISMRMRD::Acquisition*>(v_acq);
	return acq->read_dir()[idx];
}

float
jl_ismrmrd_header_phase_dir(void* v_acq, int idx) {
	if(!v_acq){
		return 0.;
	}
	ISMRMRD::Acquisition* acq = reinterpret_cast<ISMRMRD::Acquisition*>(v_acq);
	return acq->phase_dir()[idx];
}

float
jl_ismrmrd_header_slice_dir(void* v_acq, int idx) {
	if(!v_acq){
		return 0.;
	}
	ISMRMRD::Acquisition* acq = reinterpret_cast<ISMRMRD::Acquisition*>(v_acq);
	return acq->slice_dir()[idx];
}

float
jl_ismrmrd_header_patient_table_position(void* v_acq, int idx) {
	if(!v_acq){
		return 0.;
	}
	ISMRMRD::Acquisition* acq = reinterpret_cast<ISMRMRD::Acquisition*>(v_acq);
	return acq->patient_table_position()[idx];
}

// void*
// jl_ismrmrd_header_idx(void* v_acq) {
// 	if(!v_acq){
// 		return NULL;
// 	}
// 	ISMRMRD::Acquisition* acq = reinterpret_cast<ISMRMRD::Acquisition*>(v_acq);
// 	return reinterpret_cast<ISMRMRD::ISMRMRD_EncodingCounters*>(&acq->idx());
// }

int
jl_ismrmrd_header_user_int(void* v_acq, int idx) {
	if(!v_acq){
		return 0;
	}
	ISMRMRD::Acquisition* acq = reinterpret_cast<ISMRMRD::Acquisition*>(v_acq);
	return acq->user_int()[idx];
}

float
jl_ismrmrd_header_user_float(void* v_acq, int idx) {
	ISMRMRD::Acquisition* acq = reinterpret_cast<ISMRMRD::Acquisition*>(v_acq);
	return acq->user_float()[idx];
}

uint16_t
jl_ismrmrd_get_ec_kspace_encode_step_1(void* v_acq) {
	ISMRMRD::Acquisition* acq = reinterpret_cast<ISMRMRD::Acquisition*>(v_acq);
	return acq->idx().kspace_encode_step_1;
}
uint16_t
jl_ismrmrd_get_ec_kspace_encode_step_2(void *v_acq){
	ISMRMRD::Acquisition* acq = reinterpret_cast<ISMRMRD::Acquisition*>(v_acq);
	return acq->idx().kspace_encode_step_2;
}
uint16_t
jl_ismrmrd_get_ec_average(void* v_acq){
	ISMRMRD::Acquisition* acq = reinterpret_cast<ISMRMRD::Acquisition*>(v_acq);
	return acq->idx().average;
}
uint16_t
jl_ismrmrd_get_ec_slice(void* v_acq){
	ISMRMRD::Acquisition* acq = reinterpret_cast<ISMRMRD::Acquisition*>(v_acq);
	return acq->idx().slice;
}
uint16_t
jl_ismrmrd_get_ec_contrast(void* v_acq){
	ISMRMRD::Acquisition* acq = reinterpret_cast<ISMRMRD::Acquisition*>(v_acq);
	return acq->idx().contrast;
}
uint16_t
jl_ismrmrd_get_ec_phase(void* v_acq){
	ISMRMRD::Acquisition* acq = reinterpret_cast<ISMRMRD::Acquisition*>(v_acq);
	return acq->idx().phase;
}
uint16_t
jl_ismrmrd_get_ec_repetition(void* v_acq){
	ISMRMRD::Acquisition* acq = reinterpret_cast<ISMRMRD::Acquisition*>(v_acq);
	return acq->idx().repetition;
}
uint16_t
jl_ismrmrd_get_ec_set(void* v_acq){
	ISMRMRD::Acquisition* acq = reinterpret_cast<ISMRMRD::Acquisition*>(v_acq);
	return acq->idx().set;
}
uint16_t
jl_ismrmrd_get_ec_segment(void* v_acq){
	ISMRMRD::Acquisition* acq = reinterpret_cast<ISMRMRD::Acquisition*>(v_acq);
	return acq->idx().segment;
}
uint16_t
jl_ismrmrd_get_ec_user(void* v_acq, int idx){
	ISMRMRD::Acquisition* acq = reinterpret_cast<ISMRMRD::Acquisition*>(v_acq);
	return acq->idx().user[idx];
}

void
jl_ismrmrd_resize(void* v_acq, uint16_t num_samples,
	uint16_t active_channels, uint16_t traj_dimensions) {
    if(!v_acq)
	return;

    ISMRMRD::Acquisition* acq =
	reinterpret_cast<ISMRMRD::Acquisition*>(v_acq);

    acq->resize(num_samples,active_channels,traj_dimensions);
    return;
}

size_t
jl_ismrmrd_get_number_of_data_elements(void* v_acq) {
    if(!v_acq)
	return 0;

    ISMRMRD::Acquisition* acq =
	reinterpret_cast<ISMRMRD::Acquisition*>(v_acq);
    return acq->getNumberOfDataElements();
}

size_t
jl_ismrmrd_get_number_of_traj_elements(void* v_acq) {
    if(!v_acq)
	return 0;

    ISMRMRD::Acquisition* acq =
	reinterpret_cast<ISMRMRD::Acquisition*>(v_acq);
    return acq->getNumberOfTrajElements();
}

size_t
jl_ismrmrd_get_data_size(void * v_acq) {
    if(!v_acq)
	return 0;

    ISMRMRD::Acquisition* acq =
	reinterpret_cast<ISMRMRD::Acquisition*>(v_acq);
    return acq->getDataSize();
}

size_t jl_ismrmrd_get_traj_size(void * v_acq) {
    if(!v_acq)
	return 0;

    ISMRMRD::Acquisition* acq =
	reinterpret_cast<ISMRMRD::Acquisition*>(v_acq);
    return acq->getTrajSize();
}

/*
ISMRMRD_AcquisitionHeader&
jl_ismrmrd_get_head(void* v_acq) {

    ISMRMRD::Acquisition* acq =
	reinterpret_cast<ISMRMRD::Acquisition*>(v_acq);

    return acq->getHead();
}

void
jl_ismrmrd_set_head(void* v_acq, ISMRMRD_AcquisitionHeader& header) {
    if(!v_acq)
	return;

    ISMRMRD::Acquisition* acq =
	reinterpret_cast<ISMRMRD::Acquisition*>(v_acq);
    acq->setHead(header);
    return;
}
*/

complex_float_t*
jl_ismrmrd_get_data_pointer(void * v_acq) {

    ISMRMRD::Acquisition* acq =
	reinterpret_cast<ISMRMRD::Acquisition*>(v_acq);

    return acq->getDataPtr();
}

float*
jl_ismrmrd_get_traj_pointer(void * v_acq) {

    ISMRMRD::Acquisition* acq =
	reinterpret_cast<ISMRMRD::Acquisition*>(v_acq);

    return acq->getTrajPtr();
}


int
jl_ismrmrd_is_flag_set(void* v_acq, const uint16_t val) {
    ISMRMRD::Acquisition* acq =
	reinterpret_cast<ISMRMRD::Acquisition*>(v_acq);

    int b = acq->isFlagSet(val) ? 1 : 0;

    return b;
}

void
jl_ismrmrd_set_flag(void* v_acq, const uint16_t val) {
    ISMRMRD::Acquisition* acq =
	reinterpret_cast<ISMRMRD::Acquisition*>(v_acq);
    acq->setFlag(val);
    return;
}

void
jl_ismrmrd_clear_flag(void* v_acq, const uint16_t val) {
    ISMRMRD::Acquisition* acq =
	reinterpret_cast<ISMRMRD::Acquisition*>(v_acq);
    acq->clearFlag(val);
    return;
}

void jl_ismrmrd_clear_all_flags(void* v_acq) {
    ISMRMRD::Acquisition* acq =
	reinterpret_cast<ISMRMRD::Acquisition*>(v_acq);
    acq->clearAllFlags();
    return;
}

int
jl_ismrmrd_is_channel_active(void* v_acq, uint16_t channel_id) {
    ISMRMRD::Acquisition* acq =
	reinterpret_cast<ISMRMRD::Acquisition*>(v_acq);
    int b = acq->isChannelActive(channel_id) ? 1 : 0;
    return b;
}

void jl_ismrmrd_set_channel_active(void* v_acq, uint16_t channel_id) {
    ISMRMRD::Acquisition* acq =
	reinterpret_cast<ISMRMRD::Acquisition*>(v_acq);
    acq->setChannelActive(channel_id);
    return;
}

void jl_ismrmrd_set_channel_not_active(void* v_acq, uint16_t channel_id) {
    ISMRMRD::Acquisition* acq =
	reinterpret_cast<ISMRMRD::Acquisition*>(v_acq);
    acq->setChannelNotActive(channel_id);
    return;
}

void jl_ismrmrd_set_all_channels_not_active(void* v_acq) {
    ISMRMRD::Acquisition* acq =
	reinterpret_cast<ISMRMRD::Acquisition*>(v_acq);
    acq->setAllChannelsNotActive();
    return;
}
