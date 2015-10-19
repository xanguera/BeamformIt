#ifndef TDOA_H
#define TDOA_H

#include <vector>
#include "support.h"
#include "sndfile.h"
#include "fileinout.h"
//#include "delaysum.h"
//#include "boost/program_options.hpp"

#include "ffft/FFTReal.h"

//namespace po = boost::program_options;
using namespace std;

class TDOA
{
public:

	TDOA(FileInOut* fileInOut, Configuration* config);
	void init(int numChannels, float UEMGap);
	void Optimum_delays_acoustic_modelling();
	void Optimum_delays_filtering();
	void continuity_filter();
	void compute_channels_xcorr(long counter, float** chanData);
	void print_delays(int sampleRateInMs);
	float find_nbest_maximums(int *delays, float *values, int amount_max, float *xcorr_value, int margin, int fftwindow_size);

    float xcorrelation_FFTReal(int *delays, float *values, int amount_max, float *chan_data, float *ref_data, int margin, int window, int fftwindow);
    vector<int> xcorrelation_FFTReal_full(long counter, vector<vector<float> > chanData);

    void Optimum_delays_write_out();
	vector<int> get_bestChDelays(int numDelay);

	//setters and getters
	void set_numCh(int numCh){m_numCh = numCh;};
	int get_totalNumDelays(){return(m_totalNumDelays);};
	void set_totalNumDelays(int totalNumDelays){m_totalNumDelays = totalNumDelays;};
	int get_finalDelays(int channel, int delNum, int nBest){return(m_finalDelays[channel][delNum][nBest]);};
	int get_delayFilters(int channel, int frame){return(m_delayFilters[channel][frame]);};
	void set_UEMGap(float UEMGap){m_UEMGap = UEMGap;};

private:
	FileInOut *m_fileInOut;
	Configuration *m_config;
	LogInfo m_log;
	
	///keeps track of the filtering applied to each delay
	std::vector<std::vector<int> > m_delayFilters;
	std::vector<int> m_globalDelayFilters; ///same, but affecting all dimensions (is it used?)

	int m_numCh;
	int m_totalNumDelays; ///total number of delays to be computed over all the process 
	int m_fftWindowFrames; ///window used for the fft. It is the inmediate power of 2 bigger than m_windowSize
	float m_UEMGap;

	///keeps track of the best delays in the N-best for each microphone
	std::vector<std::vector<unsigned int> > m_bestClass, m_bestClass2; //indicates which of the nbest values is chosen in Viterbi
	std::vector<std::vector<std::vector<int> > > m_chDelays; 	///computed delays for channels-> positions-> N-best
	std::vector<std::vector<std::vector<float> > > m_chXcorrValues;	///computed delays for channels-> positions-> N-best
	///to hold the selected optimum delays and xcorr values
	std::vector<std::vector<std::vector<int> > > m_finalDelays;
	std::vector<std::vector<std::vector<float> > > m_finalXcorrValues; ///holds the different xcorrelation values in the maximum points



	float determine_noise_threshold();
	void noise_filter(float threshold);
	void Nbest_IndivCh_decoding();
	void Global_delays_decoding();
	void Compute_overlap_models();
	void segment_Nbestdelays();
	void segment_Nbestdelays_2best();
	void segment_Nbestdelays_multichannel(int Nbest, int block_size, int offset);
	std::vector<vector<double> > compute_trans_prob(int frame, int channel_count);
	std::vector<vector<double> > compute_trans_prob_multichannel(int frame, vector<vector<int> > combinations);
	std::vector<double> compute_xcorr_prob(int frame, int channel_count);
	std::vector<double> compute_xcorr_prob_multichannel(int frame, vector<vector<int> > combinations);
	std::vector<vector<int> > create_delays_positions_block(int Nbest, int block_size, int offset, int frame);

};










#endif
