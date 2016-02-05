#ifndef GLOBAL_H
#define GLOBAL_H

#include <vector>
#include <string>

using namespace std;

//filtering actions
#define F_NOISE 1 // applied noise filter
#define F_CONTIN 2 //applied continuity filter
#define F_OVERLAP 4 //overlap detected on this frame
#define F_UNPURE 8 //unpure segment

//variables used in the Viterbi decoding module
#define LAST 0
#define CURR 1

//maximum amount of channels	
#define MAXNUMCH 128

//configuration parameters relative to the delay&sum algorithm, not the data
// Parameters for the data are in the DS class
struct Configuration
{
    string m_config_file; //where is the configuration file located
    int delay_variance;
    int ovl_variance;
    int INDIV_CHANNELS;
    int OUT_FORMAT;
    int COMPUTE_REFERENCE;
    int DO_OPTIMUM_DELAYS;
    int DO_AVOID_BAD_FRAMES;
    int badFramesRatio; //parameter in the threshold to determine a bad frame from a good one, the higher the better
    int DO_USE_UEM_FILE;
    int DO_ADAPT_WEIGHTS;
    int DO_WRITE_SPH_FILES;
    int DO_ACOUSTIC_MODELLING;
    int DO_NOISE_THRESHOLD;
    int DO_OUTPUT_RESIDUAL; //flag to output the residual between 2 signals
    int DO_SIGNAL_QUALITY; //flag to compute a signal quality measure from the output file
    int DO_COMPUTE_SKEW; //flag to whether we compute the skew between signals
    int DO_OUTPUT_OVERLAP; //flag whether to compute possile overlaps from the audio
    int DO_DELAYS_PADDING; //flag whether to add extra delays to better cover all acoustics

    string AUDIOROOT; //where to look for the audio files
    string RESULTPATH; //where we write the resulting files
    string CHANNELSFILE; //were to find the config file
    string UEMFILE; //where to find the UEM file
    string SHOWNAME; //show name

    int rate;
    int window; //analysis window, in ms
    int windowFrames; //analysis window, in frames

    float min_xcorr; // minimum xcorr value to consider it not noise

    int margin; //area in ms around where we seach for the max of the autocorrelation
    int marginFrames; //same as margin, but in number of frames
    int extra_delays_padding; //number of repeated delays appended to the end in the output file

    int skewMargin; //area in ms around where we search for the skew between signals

    int m_printLevel; //print level for the log
  unsigned int nbest_amount; //amount of Nbest delays considered
  int reference_channel;

  int print_features;

  //acoustic modelling parameters
  float trans_weight_nbest;
  float trans_weight_multi;

  float noise_percent;

  int amount_model_feats;

  //int full_path;

  int output_second_best_sph;
  
	float refAdaptRatio; //reference channel selection adaptation ratio
	
};


#endif
