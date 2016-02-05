#ifndef FILEINOUT_H
#define FILEINOUT_H

#include "global.h"
#include "support.h"
#include "sndfile.h"

#include <string>

class FileInOut
{

public:

  //input files
  SNDFILE* inputfd[MAXNUMCH];
  SF_INFO file_info;
  SF_INFO file_input_info; //saves the original file input information
  char input_file_extension[4];
  //output files
  SNDFILE* outfd[2]; //1st and 2nd best
  SNDFILE* outputfd[MAXNUMCH];
  //delays
  FILE* delfd;
  FILE* delfd2;  
  FILE* ovlfd; //overlap
  FILE* infofd; //general info file
  FILE* featwfd; //weight feature file

	FileInOut(Configuration* config);
	~FileInOut();
	void init();
	int readChanData(int channel, float* chanData, long startSample, int numSamples);
	void delays_to_file(vector<vector<vector<int> > > & finalDelays, vector<vector<vector<float> > > & finalXcorrValues, float UEMGap);

	void Open_Output_Channels();
	int Open_Input_Channels();
	void Close_Channels();
	void Channels_Sum(int best_out);

	//getters and setters
	int get_numCh(){return(m_numCh);}
	long get_nFrames(){return(m_frames);};
	long get_nFramesInit(){return(m_framesInit);};
	long get_sampleRate(){return(m_sampleRate);};
	int get_sampleRateInMs(){return(m_sampleRateInMs);};
	int get_framesOut(){return(m_framesOut);};
	char* get_fileOut(){return(m_fileOut);};

private:

	char* m_fileOut; //name of the beamformed file for the 1st best delays

	LogInfo m_log;	///log and Info class
	Configuration *m_config;

	long m_frames; ///number of total frames in the region of the files being processed
	long m_framesInit; ///number of frames in total in the input files
	long m_sampleRate; ///sample rate of the files
	int m_sampleRateInMs; /// samplerate in milliseconds
	///endings of all channels to do the delay and sum with
	vector<char *> m_channels; //path to each of the audio input channels
	///number of channels/wavefiles
	int m_numCh;	 
	///number of files containing such channels (some files might contain multiple channels
	int m_numFiles;
	int m_framesOut; ///frames to write out when doing the sum process

};




#endif
