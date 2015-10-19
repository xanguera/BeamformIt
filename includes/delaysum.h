#ifndef DELAYSUM_H
#define DELAYSUM_H

#include "global.h"
#include "fileinout.h"
#include "support.h"
#include "tdoa.h"
//#include "boost/program_options.hpp"

//namespace po = boost::program_options;
using namespace std;

class DelaySum {

public:

	//variables regarding the channels I have to process
	int m_numCh; ///number of channels available for the show

	int totalNumDelays; //number of TDOA values being computed

	int m_windowSize; ///size in frames

	vector<float> m_refData;	///pointer to the reference channel data
	vector< vector<float> > m_chanData;	///pointer to the channels data

	int m_framesOut; ///frames to write out when doing the sum process

	float m_overallWeight;	///overall weigthing factor
	float m_channelWeight[MAXNUMCH];	///temp variables for the sum of the channels
	float m_outWeightSNR[MAXNUMCH]; ///for the snr modification (unused??)
	float m_globalWeight[MAXNUMCH]; ///adaptable global weights values

	///holds the pairwise xcorr values for a given frame, with a given optimum delay
	double m_localXcorr[MAXNUMCH][MAXNUMCH]; ///holds the local xcorrelation of delayed channels
	int m_chanFramesRead[MAXNUMCH];///frames read from the file for a given frame
	double m_localEnergy[MAXNUMCH];///holds the local energy value computed for each channel 

	long m_frames; ///number of total frames in the region of the files being processed
	long m_framesInit; ///number of frames in total in the input files
	long m_sampleRate; ///sample rate of the files
	int m_sampleRateInMs; /// samplerate in milliseconds
	float m_UEMGap; ///UEM offset start from 0
	float m_UEMEndTime; ///UEM end time
    vector<int> m_skew; ///skew for each of the channels
	int m_biggestSkew; ///biggest of the skews
	float m_ave_val_segment[MAXNUMCH]; ///maximum average value for each signal 

	char m_dataSet[64]; ///char string containing the site's origin of the recording (i.e. ICSI, NIST...)

	int m_printLevel; ///sets the level(inclusive) up to where print information should be displayed, the higher the more strict

	~DelaySum();
	DelaySum(Configuration *config, TDOA *tdoa, FileInOut* fileInOut);
	void init();
	void get_channels_data(long file_start);
 	void get_channels_data(long file_start, vector<int> ch_delays);
	void GetChInfo();
    void Set_Show_Parameters();
    //void input_parameters(int argc, char *argv[], po::variables_map &vm);
	void compute_scalling_factor();
	void compute_skew();
	void SetUEMlimits();
	void Channels_Sum(int best_out);
	void Set_Channels_Weights();
	long get_start_frame(long counter);
	int find_best_channel_xcorr();
	int find_best_channel_xcorr_adapt(long start_frame, vector<int> ch_delays);
	int find_best_channel_xcorr_adapt();
	void computeChannelsDifference();

	//getters and setters
	long get_frames(){return(m_frames);};
	long get_sampleRate(){return(m_sampleRate);};
	int get_windowSize(){return(m_windowSize);};
	int get_skew(int channel){return(m_skew[channel]);};
	vector<float> get_chanData(int channel){return(m_chanData[channel]);};
	float get_UEMGap(){return(m_UEMGap);};
	void set_numCh(int numCh){m_numCh = numCh;};
	void set_UEMGap(float value){m_UEMGap = value;};
	void signalQuality();
	void computeTDOAValues(long start_frame, int delayNum);

private:

	LogInfo m_log;	///log and Info class
	TDOA* m_tdoa; //tdoa data and processing functions
	Configuration* m_config; //config structure
	FileInOut* m_fileInOut; //input output class

	int skew_detector(SNDFILE* file_in, SNDFILE* file_ref, long frames, long samplerate);
	float delay_detector(SNDFILE* file_in1, SNDFILE* file_in2, long frames, long samplerate, int numNbest);
	float variance_detector(SNDFILE* file_in, SNDFILE* file_ref, long frames, long samplerate, int amount_pieces, int do_filter);
    void GetUEMInfo(char *uemFile, char *showName, float *start_time, float *end_time);
	//void delays_to_file(FileInOut &fileInOut, Configuration &config);


	//functions for channels sum
	void Compute_Sum_Weights(int frame);
	void compute_local_xcorr_values(vector<vector<float> > &chan_out);
	void compute_local_energy(vector<vector<float> > &chan_out);


};


#endif
