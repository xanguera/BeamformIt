/* delaysum.cpp, Implements the delay_sum and other related techniques */
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <vector>
#include <string>
#include <set>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <cstring>

//#include "samplerate.h"
#include "sndfile.h"
#include "global.h"
#include "delaysum.h"
#include "support.h"

//namespace po = boost::program_options;

//char tmp_string[1024]; //used to output messages

/*!
  Default destructor
*/
DelaySum::~DelaySum()
{
}

/*!
  Default constructor
*/
DelaySum::DelaySum(Configuration *config, TDOA *tdoa, FileInOut* fileInOut)
{
	//assign the data structs
	m_tdoa = tdoa;
	m_config = config;
	m_fileInOut = fileInOut;

	//initialize variables to default values
	m_biggestSkew = 0;
    m_skew.clear();
    m_skew.resize(MAXNUMCH, 0);
}

/*!
  Initialize all necessary variables in the delay_sum class, before using them.
  \return Nothing
*/

void DelaySum::init()
{

  //define the window size in frames
  (*m_config).windowFrames = (int)((*m_config).window *(m_sampleRateInMs)); //the window we compute in frames

  //we reserve the memory for the data
	m_chanData.resize(m_numCh);
  for(int count=0; count< m_numCh; count++)
  {
    m_chanData[count].resize((*m_config).windowFrames, 0);
  }

  m_refData.resize((*m_config).windowFrames, 0);

  //we compute how many delays we will process
  //NOTE: we are taking a flooring of the resulting value. This means we will be missing a bit of the end of the audio,
  //      which should not be dramatic for most applications
  int totalNumDelays = (int)((m_frames - (*m_config).windowFrames - m_biggestSkew - m_UEMGap)/((*m_config).rate*m_sampleRateInMs));
  (*m_tdoa).set_totalNumDelays(totalNumDelays);

  char tmp_string[1024];
	sprintf(tmp_string, "Total number of delays to be computed: %d\n", totalNumDelays);
  m_log.print_this(tmp_string, 1);
}

/*!
   Used to process all input parameters and fill the appropriate config values with default and input values
   \param argc Is the number of elements in the command line
   \param *argv[] Is the command line to be processed
   \param &vm Variables_map used by the library program_options to hold all input parameters read from the command line or the config file. All these parameters are transfered to the config structure afterwards.
   \return Nothing
*/

/*!
  defines some show-dependent parameters
  \param &vm Variables_map used by the library program_options to hold all input parameters read from the command line or the config file. All these parameters are transfered to the config structure afterwards.
  \param FileInOut &fileInOut class pointer with input/output file descriptors
  \return Nothing
 */
void DelaySum::Set_Show_Parameters()
{
    //set log level
    m_log.set_printLevel((*m_config).m_printLevel);

    /*
    //perform some checks on parameters
    //check the use of INDIV_CHANNELS and full_path (are incompatible)
    if((*m_config).full_path && (*m_config).INDIV_CHANNELS)
    {
      if((*m_config).INDIV_CHANNELS == 1)
      {
        m_log.print_this("WARNING: Use of individual channel output when using full path input files is not supported\n",5);
        m_log.print_this("  Changing to one channel only output\n",5);
        (*m_config).INDIV_CHANNELS = 0;
      }
    }
    */
  //Set the channels necessary information
  m_frames = (*m_fileInOut).get_nFrames();
  m_framesInit = (*m_fileInOut).get_nFramesInit();
  m_sampleRate = (*m_fileInOut).get_sampleRate();
  m_sampleRateInMs = (*m_fileInOut).get_sampleRateInMs();
  m_numCh = (*m_fileInOut).get_numCh();
  m_framesOut = (*m_fileInOut).get_framesOut();


  //define the margins for the input speech depending one the data set
  //float margin; //in ms, used +-margin

  /* NOTE: now I just consider a single value
  //get the show name
  int count=0;
  int tmp_margin;

  while((*m_config).SHOWNAME[count] != '_' && (*m_config).SHOWNAME[count] != '\0')
  {
    m_dataSet[count] = (*m_config).SHOWNAME[count];
    count++;
  }
  m_dataSet[count]='\0';

  if(!strcmp(m_dataSet,"ICSI"))
  {
    tmp_margin = 20; // ~ 12(skew) + 8(mics) 20
  }
  else if(!strcmp(m_dataSet,"NIST"))
  {
    tmp_margin = 10; // ~ 10(mics)
  }
  else if(!strcmp(m_dataSet,"CMU"))
  {
    tmp_margin = 15;
  }
  else if(!strcmp(m_dataSet,"LDC"))
  {
    tmp_margin = 10;
  }
  else
  {
    char tmp_string[1024]; sprintf(tmp_string,"Unknown data set: %s , setting default values\n",m_dataSet);
    m_log.print_this(tmp_string, 1);
    tmp_margin = 20;
  }

  //we set the values if they are not set by configuration
  if((*m_config).margin < 0)
  {
    (*m_config).margin = tmp_margin;
  }
*/

  //convert the margin in ms to number of frames
  (*m_config).marginFrames = (int)((float)((*m_config).margin * m_sampleRateInMs));
  char tmp_string[1024]; sprintf(tmp_string,"Marginst for delays in frames: %d in ms: %d\n",(*m_config).marginFrames, (*m_config).margin);
  m_log.print_this(tmp_string, 1);

}


/*!
  Compute the average maximum xcorr between all channel pairs and determine which one is best
  \return The best/reference channel according to xcorr measures.
*/

int DelaySum::find_best_channel_xcorr()
{
  vector<vector<float> > xcorr_values;
	xcorr_values.resize(m_numCh);
  for(int i=0; i<m_numCh; i++)
		xcorr_values[i].resize(m_numCh, 0);

  //first we compute all pair-wise channels average xcorr
  m_log.print_this("Determining the reference channel via average xcorr values\n", 1);
  for(int count1=0; count1< m_numCh-1; count1++)
  {
    for(int count2=count1+1; count2<m_numCh; count2++)
    {
      char tmp_string[1024]; sprintf(tmp_string,"Computing the average xcorr between %d and %d\n", count1, count2);
      m_log.print_this(tmp_string, 1);
      xcorr_values[count1][count2] = delay_detector((*m_fileInOut).inputfd[count1], (*m_fileInOut).inputfd[count2], m_framesInit, m_sampleRate, 2);
      xcorr_values[count2][count1] = xcorr_values[count1][count2];
    }
  }

  //now we accumulate all xcorrs for each channel and determine the biggest one
  vector<float> ave_xcorr(m_numCh, 0);
  for(int count1=0; count1< m_numCh; count1++)
  {
    for(int count2=0; count2<m_numCh; count2++)
    {
      ave_xcorr[count1] += xcorr_values[count1][count2];
    }
    ave_xcorr[count1] /= m_numCh;
    char tmp_string[1024]; sprintf(tmp_string,"Average xcorr for channel %d: %f\n", count1, ave_xcorr[count1]);
    m_log.print_this(tmp_string, 0);
  }

  int best_channel = 0;
  float best_xcorr = 0;
  for(int i=0; i<m_numCh; i++)
  {
    if(ave_xcorr[i] > best_xcorr)
    {
      best_xcorr = ave_xcorr[i];
      best_channel = i;
    }
  }
  char tmp_string[1024]; sprintf(tmp_string,"Reference channel set to channel %d\n", best_channel);
  m_log.print_this(tmp_string, 0);
  return(best_channel);
}

/*!
	Compute a measure of the quality of the beamformed signal by performing the average of the xcorr values of such signal
	with all individual input signals. Such metric is not comparable across recordings, but is useful to compare two different
	beamforming outputs for one particular recording. The higher the value, the better
*/
void DelaySum::signalQuality()
{
  vector<float> xcorr_values(m_numCh, 0);

	//create a pointer to the output file for reading
	SNDFILE* BF;
	SF_INFO file_info;
	file_info.format = 0; //mandatory when opening for reading
	char* fileOut = (*m_fileInOut).get_fileOut(); //get the name of the file
	BF = sf_open(fileOut, SFM_READ, &(file_info));

  //first we compute the average xcorr for all channels
  m_log.print_this("Computing the xcorr for all channels and the beamformed signal\n", 1);
	fprintf((*m_fileInOut).infofd, "Quality score per channel: ");
  for(int count1=0; count1< m_numCh; count1++)
	{
    xcorr_values[count1] = delay_detector(BF, (*m_fileInOut).inputfd[count1], m_framesInit, m_sampleRate, 1);
		fprintf((*m_fileInOut).infofd, "%f ", xcorr_values[count1]);
  }
	fprintf((*m_fileInOut).infofd, "\n");

	//make the average on these averages
	float ave_xcorr = 0;
	for(int count1=0; count1< m_numCh; count1++)
		ave_xcorr += xcorr_values[count1];

	ave_xcorr /= m_numCh;

  char tmp_string[1024]; sprintf(tmp_string,"Quality score for beamformed signal is %f\n", ave_xcorr);
  m_log.print_this(tmp_string, 0);

	//write out the score into the info file
	fprintf((*m_fileInOut).infofd, "Overall quality score: %f\n", ave_xcorr);
}



/*
  Calls find_best_channel_xcorr for the initial case
  \return the best channel number
*/
int DelaySum::find_best_channel_xcorr_adapt()
{
  vector<int> tmp_delays(m_numCh, 0);
  return(find_best_channel_xcorr_adapt(0, tmp_delays));
}

/*
  Version of the algorithm above to get the best channel considering a single window for each channel
	for the xcorr comparison but relies in an adaptative channel weighting to select the optimum channel
  \param dtart_frame starting frame for the computation
  \param ch_delays the channel optimum delays at this point
  \return The best channel found at this step
*/
int DelaySum::find_best_channel_xcorr_adapt(long start_frame, vector<int> ch_delays)
{
  //variable to hold the xcorr matrix
  vector<vector<float> > xcorr_values;//[m_numCh][m_numCh];
  for(int i=0; i<m_numCh; i++)
  {
    vector<float> tmp_vector(m_numCh, 0);
    xcorr_values.push_back(tmp_vector);
  }

  //get the channels data with the given delays
  get_channels_data(start_frame, ch_delays);

  //compute the energy of each one of the channels, stores them in m_localEnergy[]
  compute_local_energy(m_chanData);


	//Alternative: just use the energy in every channel to determine the best channel in every step
	float alpha = (*m_config).refAdaptRatio;
	int best_channel = 0;
	float best_energy = -1;
    for(int channel_count=0; channel_count< m_numCh; channel_count++)
    {
        float normEnergy = m_localEnergy[channel_count]/m_ave_val_segment[channel_count];
        m_globalWeight[channel_count] = (1-alpha) * m_globalWeight[channel_count] + alpha * normEnergy;

        //look for the best
        if(m_globalWeight[channel_count] > best_energy)
        {
            best_energy = m_globalWeight[channel_count];
            best_channel = channel_count;
        }
    }


    /*
	//Alternative 2: normalize the energy in every channel by a slow moving average
	float alpha = 0.1;
	int best_channel = 0;
	float best_energy = -1;
  for(int channel_count=0; channel_count< m_numCh; channel_count++)
  {
		float normEnergy = m_localEnergy[channel_count]/m_globalWeight[channel_count];
    m_globalWeight[channel_count] = (1-alpha) * m_globalWeight[channel_count] + alpha * m_localEnergy[channel_count];

		//look for the best
    if(normEnergy > best_energy)
    {
      best_energy = normEnergy;
      best_channel = channel_count;
    }
  }
*/

/*
	//Alternative 3: select only among the channel pair with highest xcorr between reference and other channels. Selection is done by energy
	float alpha = 0.1;
	int best_channel = 0;
	float best_energy = -1;
  for(int channel_count=0; channel_count< m_numCh; channel_count++)
  {
		float normEnergy = m_localEnergy[channel_count]/m_globalWeight[channel_count];
    m_globalWeight[channel_count] = (1-alpha) * m_globalWeight[channel_count] + alpha * m_localEnergy[channel_count];

		//look for the best
    if(normEnergy > best_energy)
    {
      best_energy = normEnergy;
      best_channel = channel_count;
    }
  }
*/

/*
  //compute the xcorr values for all combinations, stores them in m_localXcorr[]
  compute_local_xcorr_values(m_chanData);

  //from the computed local xcorr values compute the relative weights, adapt the overall weights and decide on the best channel at this point

  vector<double> vector_xcorr(m_numCh,0);
  double sum_vector_xcorr = 0;

  //sum the xcorr values of each channel with all others (except itself)
  for(int channel_count1=0; channel_count1< m_numCh; channel_count1++)
  {
    vector_xcorr[channel_count1] = 0;
    for(int channel_count2=0; channel_count2< m_numCh; channel_count2++)
    {
      //do not count itself
      if(channel_count1 != channel_count2)
      {
        vector_xcorr[channel_count1] += m_localXcorr[channel_count1][channel_count2];
      }
    }
    sum_vector_xcorr += vector_xcorr[channel_count1];
  }

  //in weird cases the sum could be 0, we avoid that (This would happen when all xcorr values are negative)
  if(sum_vector_xcorr == 0)
	{
		sum_vector_xcorr = 1;
		for(int i=0; i<m_numCh; i++)
			vector_xcorr[i] = 1.0/(double)(m_numCh);
	}

  //convert local xcorrs into weights and adapt the global weight
  double local_xcorr_weight;

	float alpha = 0.05;
  for(int channel_count=0; channel_count< m_numCh; channel_count++)
  {
    local_xcorr_weight = vector_xcorr[channel_count]/sum_vector_xcorr;
    m_globalWeight[channel_count] = (1-alpha) * m_globalWeight[channel_count] + alpha * local_xcorr_weight;
  }

  //find the best channel looking at the global weights
	//set as default value the current reference channel, in case all others have the same value as it, we just keep the same
  int best_channel = (*m_config).reference_channel;
  float best_xcorr = m_globalWeight[(*m_config).reference_channel];
  for(int i=0; i<m_numCh; i++)
  {
    if(m_globalWeight[i] > best_xcorr)
    {
      best_xcorr = m_globalWeight[i];
      best_channel = i;
    }
  }
*/

	return(best_channel);
}



/*!
  Computes the scalling factor in order to adapt the average dynamic range of the input data to the maximum at the ouput
  \return Nothing
*/
void DelaySum::compute_scalling_factor()
{
  float max_val_segment[m_numCh];
  vector<float> median_val_segment;
  int num_segments=0;
  int total_num_segments_initial = 10; //we will only compute (at max) these number of segments
  int total_num_segments = total_num_segments_initial; // it can be variable
  int signal_scroll; //advancement in the signal for the next segment
  int segment_duration = 10 * m_sampleRate; //10 seconds segments
  float sum_weighting;
  float *signal = (float *)malloc(segment_duration * sizeof(float));//10 seconds
  long start_place=0;
  long amount_frames_read=0;
  long amount_frames_total = 0;
  m_log.print_this("Calculating the maximum values of the signal, to do a signal adjustment\n", 1);
  for(int channel_count=0; channel_count< m_numCh; channel_count++)
  {
    median_val_segment.clear();
    amount_frames_total = 0;
    //max_val[channel_count] = 0;
    //rms_value = 0;
    m_ave_val_segment[channel_count] = 0;
    num_segments = 0;
    start_place = (long)(m_UEMGap);
    //check if we have enough signal for the number of segments
    if(((m_frames - start_place)/segment_duration) < total_num_segments_initial)
    {
      total_num_segments = int((m_frames - start_place) / segment_duration);
      if (total_num_segments < 1)
      {
				total_num_segments = 1;
				segment_duration = m_frames - start_place; //we do not have more frames
      }
      char tmp_string[1024]; sprintf(tmp_string,"Set the total number of segments from %d to %d for computing scalling factor, with segment duration %d\n", total_num_segments_initial, total_num_segments, segment_duration);
      m_log.print_this(tmp_string, 1);
    }
    signal_scroll = int((m_frames - start_place)/total_num_segments);
    char tmp_string[1024]; sprintf(tmp_string,"Processing channel %d\n", channel_count);
    m_log.print_this(tmp_string, 1);
    while(start_place < m_frames)
    {
      sf_seek((*m_fileInOut).inputfd[channel_count], start_place, SEEK_SET);
      amount_frames_read = sf_readf_float((*m_fileInOut).inputfd[channel_count], signal, segment_duration);//10 seconds
      amount_frames_total += amount_frames_read;
      max_val_segment[channel_count] = 0;
      for(int i=0; i<amount_frames_read; i++)
          if(max_val_segment[channel_count] < fabs(signal[i]))
              max_val_segment[channel_count] = fabs(signal[i]);
      median_val_segment.push_back(max_val_segment[channel_count]);
      num_segments++;
      start_place += signal_scroll;
    }

    sort(median_val_segment.begin(), median_val_segment.end());
    m_ave_val_segment[channel_count] = median_val_segment[(int)(total_num_segments/2)];
    sprintf(tmp_string,"The Median maximum energy for channel %d is %f\n",channel_count, m_ave_val_segment[channel_count]);
    m_log.print_this(tmp_string, 1);
  }

  sum_weighting = 0;
  for(int channel_count=0; channel_count< m_numCh; channel_count++)
  {
    //we need to make it relative to the channel's weight (normally 1/m_numCh)
    sum_weighting += m_ave_val_segment[channel_count];
  }

  // now the maximum value for the input data is 1.0
  m_overallWeight = (0.3*m_numCh)/sum_weighting;

  char tmp_string[1024]; sprintf(tmp_string,"Weighting calculated to adjust the signal: %f\n",m_overallWeight);
  m_log.print_this(tmp_string, 1);
  free(signal);
}

/*!
  Computes the skew between channels in the recordings. This was initially created to
  solve an error in the ICSI meetings where an error in the recording equipment caused
  the different channels to start recording at different instants. Right now it is available
  through an input parameter and can be used to perform an initial alignment between channels that
  have not started at the same time.
  \return Nothing
*/
void DelaySum::compute_skew()
{
    m_biggestSkew = 0;
    for(int count=0; count< m_numCh; count++)
    {
        m_skew[count] = 0; //by default is 0 for all meetings
        if(count != (*m_config).reference_channel)
        {
            //I'm interested in the skew over ALL the meeting, not only on the 10 minutes that are spoken
            char tmp_string[1024]; sprintf(tmp_string,"Finding the skew for channel %d\n", count);
            m_log.print_this(tmp_string, 1);
            m_skew[count] = skew_detector((*m_fileInOut).inputfd[count], (*m_fileInOut).inputfd[(*m_config).reference_channel], m_framesInit, m_sampleRate);
            if(m_skew[count] > m_biggestSkew) { m_biggestSkew = m_skew[count];}
        }
    }
}

/*!
  Computes the skew for one channel pair (only used in the ICSI channels)
  \param file_in SNDFILE pointer to the channel input
  \param file_ref SNDFILE pointer to the reference channel
  \param frames Number of frames in the channels to be processed
  \param samplerate Sample rate of the data
  \return the skew found
 */

int DelaySum::skew_detector(SNDFILE* file_in, SNDFILE* file_ref, long frames, long samplerate)
{
  //we select 25 pieces of 10 seconds each and calculate the average delay (I hope this way I'll get different people talking)

  int number_pieces = 25;
  int count;
  long scroll = frames/(number_pieces+2);
  long length = samplerate * 20; //20 seconds
  int maxSkew = (int)((*m_config).skewMargin * samplerate/1000); //maximum margin
  //make sure the FFT length is a power of 2
  long FFTSize = int(pow(2.0, ceil(log(double(length))/log(2.0))));

  float *data_in = (float*)malloc(length*sizeof(float));
  float *data_ref = (float*)malloc(length*sizeof(float));
  int delays;
  float values;
  float average_skew = 0;

  for(count=1; count<=number_pieces; count++)
  {
    //we position the pointer in the place
    sf_seek(file_in, scroll*count, SEEK_SET);
    sf_seek(file_ref, scroll*count, SEEK_SET);
    //we obtain the samples
    if(sf_readf_float(file_in, data_in, length) != length)
    {
      m_log.print_this("ERROR: Error retrieving data for file_in, not enough samples in the file\n", 10);
      exit(1);
    }

    if(sf_readf_float(file_ref, data_ref, length) != length)
    {
      m_log.print_this("ERROR: Error retrieving data for file_ref, not enough samples in the file\n", 10);
      exit(1);
    }
    //we do the xcorr
    (*m_tdoa).xcorrelation_FFTReal(&delays, &values, 1, data_in, data_ref, maxSkew, length, FFTSize);
    average_skew += delays;
    char tmp_string[1024]; sprintf(tmp_string,"Skew calculation iteration %d: %d (%f)\n", count, delays, values);
    m_log.print_this(tmp_string, 0);
  }
  average_skew /= number_pieces;
  char tmp_string[1024]; sprintf(tmp_string,"Average skew: %f\n", average_skew);
  m_log.print_this(tmp_string, 1);

  free(data_in);
  free(data_ref);
  return int(average_skew);
}


/*!
  looks into the UEM file and retrieve the start and end time for the defined show
  \param uemFile Name of the file to look into
  \param show Name of the show/meeting being processed
  \param start_time Start time to be processed (in seconds)
  \param end_time End time to be processed (in seconds)
  \return Nothing
 */

void DelaySum::GetUEMInfo(char *uemFile, char *showName, float *start_time, float *end_time)
{

  FILE *uemfd;
  uemfd = fopen(uemFile, "r");
  if(uemfd == NULL)
  {
    char tmp_string[1024]; sprintf(tmp_string,"Error Opening file %s\n",uemFile);
    m_log.print_this(tmp_string, 10);
    exit(1);
  }

  char line[256];
  while(fgets(line, 300, uemfd) != NULL)
  {
    //check if it's the line we want
    if(!strncmp(line, showName, strlen(showName)))
    {
      sscanf(line,"%*s %*d %f %f",start_time, end_time);
    }
  }
  fclose(uemfd);
  return;
}


/*!
  Sets the time limits of processing in a file
  \param config Configuration structure holding condifuration parameters used throughout the system.
  \return Nothing
 */

void DelaySum::SetUEMlimits()
{
  m_UEMGap = 0; //in frames
  //if(strcmp((*m_config).UEMFILE,"none"))
  if((*m_config).UEMFILE != "")
  {
    //we are using a UEM file
    GetUEMInfo((char*)(*m_config).UEMFILE.c_str(), (char*)((*m_config).SHOWNAME.c_str()), &m_UEMGap, &m_UEMEndTime);
    char tmp_string[1024]; sprintf(tmp_string,"UEM Start time %f , end time %f\n", m_UEMGap, m_UEMEndTime);
    m_log.print_this(tmp_string, 1);

    if(m_frames >=(long)(m_UEMEndTime * m_sampleRate))
    {
      m_frames = (long)(m_UEMEndTime * m_sampleRate); //define when to finish
    }
    else
    {
      m_log.print_this("WARNING: number of frames in raw file are fewer than claimed in UEM file, proceeding with real #\n", 5);
    }

    m_UEMGap = m_UEMGap*m_sampleRate;
  }

  char tmp_string[1024]; sprintf(tmp_string,"Processing frames from %.0f to %ld\n",m_UEMGap, m_frames);
  m_log.print_this(tmp_string, 1);

}


/*!
  Computes the average of maximum xcorr values between 2 signals
  I.e. the average delay between two signals
  \param file_in1 Comparison file 1
  \param file_in2 Comparison file2
  \param frames number of frames in each file
  \param samplerate Sample Rate of the data
	\param numNbest Number of xcorr Nbest values I consider
  \return the delay average
*/

float DelaySum::delay_detector(SNDFILE* file_in1, SNDFILE* file_in2, long frames, long samplerate, int numNbest)
{

  int amount_pieces = 200;
  int count;
  long scroll = frames/(amount_pieces+2);
  long length = samplerate * 1; //1 seconds
  //make sure the FFT length is a power of 2
  long FFTSize = int(pow(2.0, ceil(log(double(length))/log(2.0))));
  int margin = (int)((*m_config).margin * samplerate/1000); //we give ample margin to find good correlations
  float *data_in = (float*)malloc(FFTSize*sizeof(float));
  float *data_ref = (float*)malloc(FFTSize*sizeof(float));
  int delays[2];
  float values[2];
  float average_xcorr = 0;
  int amount_computed_pieces = 0;

  //char tmp_string[1024]; sprintf(tmp_string, "Amount of scroll: %ld\n", scroll);
  //m_log.print_this(tmp_string, 1);

  for(count=1; count<=amount_pieces; count++)
  {
    amount_computed_pieces++;
    //check that there is enough data
    if(count*scroll+FFTSize<frames)
    {
      //we position the pointer in the place
        //printf("Comparing from frame %d to %d\n", scroll*count, scroll*count+length);
      sf_seek(file_in1, scroll*count, SEEK_SET);
      sf_seek(file_in2, scroll*count, SEEK_SET);
      //we obtain the samples
      if(sf_readf_float(file_in1, data_in, length) != length)
      {
        m_log.print_this("ERROR: Error retrieving data for file_in1, not enough samples in the file\n", 10);
        exit(1);
      }
      if(sf_readf_float(file_in2, data_ref, length) != length)
      {
        m_log.print_this("ERROR: Error retrieving data for file_in2, not enough samples in the file\n", 10);
        exit(1);
      }
      //we do the xcorr

      (*m_tdoa).xcorrelation_FFTReal(delays, values, numNbest, data_in, data_ref, margin, length, FFTSize);

      for(int i=0; i<numNbest; i++)
          average_xcorr += values[i];
      char tmp_string[1024]; sprintf(tmp_string, "xcorr calculation iteration %d: %f\n", count, values[0]);
      m_log.print_this(tmp_string, 0);
    }
  }
  average_xcorr /= (numNbest * amount_computed_pieces);
  char tmp_string[1024]; sprintf(tmp_string, "Average xcorr final: %f\n", average_xcorr);
  m_log.print_this(tmp_string, 1);

  free(data_in);
  free(data_ref);
  return(average_xcorr);

}

/*!
  we select "amount_pieces" pieces of 5 seconds each and calculate the average delay (I need this time only one person at a time talking)
  if "do_filter" is 1, we don't count delays with xcorr< 0.1

  \param file_in Comparison file pointer
  \param file_ref Reference file pointer
  \param frames
  \param samplerate
  \param amount_pieces
  \param do_filter
  \return The variance of the signal computed

*/

float DelaySum::variance_detector(SNDFILE* file_in, SNDFILE* file_ref, long frames, long samplerate, int amount_pieces, int do_filter)
{

  int count;
  long scroll = frames/(amount_pieces+2);
  long length = samplerate * 5; //5 seconds
  //make sure the FFT length is a power of 2
  long FFTSize = int(pow(2.0, ceil(log(double(length))/log(2.0))));
  int margin = (int)((*m_config).margin * samplerate/1000); //we give ample margin to find good correlations
  float *data_in = (float*)malloc(length*sizeof(float));
  float *data_ref = (float*)malloc(length*sizeof(float));
  int delays;
  float values;

  float average_delay = 0;
  float average_delay2 = 0;
  float variance =0;

  float average_xcorr=0;

  int used_pieces = 0;

  char tmp_string[1024]; sprintf(tmp_string, "Amount of scroll: %ld\n", scroll);
  m_log.print_this(tmp_string, 1);

  for(count=1; count<=amount_pieces; count++)
  {
    //we position the pointer in the place
    sf_seek(file_in, scroll*count, SEEK_SET);
    sf_seek(file_ref, scroll*count, SEEK_SET);
    //we obtain the samples
    if(sf_readf_float(file_in, data_in, length) != length)
    {
      m_log.print_this("ERROR: Error retrieving data for file_in, not enough samples in the file\n", 10);
      exit(1);
    }
    if(sf_readf_float(file_ref, data_ref, length) != length)
    {
      m_log.print_this("ERROR: Error retrieving data for file_ref, not enough samples in the file\n", 10);
      exit(1);
    }
    //we do the xcorr
    (*m_tdoa).xcorrelation_FFTReal(&delays, &values, 1, data_in, data_ref, margin, length, FFTSize);

    if(do_filter)
    {
      if(values > 0.05)
      {
        average_delay += delays;
        average_delay2 += delays*delays;
        used_pieces++;
        average_xcorr += values;
        char tmp_string[1024]; sprintf(tmp_string, "Skew calculation iteration %d: %d (%f)\n", count, delays, values);
        m_log.print_this(tmp_string, 0);
      }
    }
    else
    {
      average_delay += delays;
      average_delay2 += delays*delays;
      used_pieces++;
      average_xcorr += values;
      char tmp_string[1024]; sprintf(tmp_string, "Skew calculation iteration %d: %d (%f)\n", count, delays, values);
      m_log.print_this(tmp_string, 0);
    }


  }
  average_delay /= used_pieces;
  variance = average_delay2/used_pieces - average_delay * average_delay;
  average_xcorr /= used_pieces;
  sprintf(tmp_string, "variance: %f xcorr: %f\n", variance, average_xcorr);
  m_log.print_this(tmp_string, 1);

  free(data_in);
  free(data_ref);
  return (variance);
}



/*!
  Computes the channels weights by adapting them from the xcorr values across channels and as a
  function of the previous values. It also cancels out any channel which is too bad.

  \param frame Segment number being computed
  \return Nothing
*/

void DelaySum::Compute_Sum_Weights(int frame)
{


  //  Normalize the local xcorr over all channels before weight adaptation

  //if we want to avoid bad frames to inferfere:
  double vector_xcorr[m_numCh];
  double sum_vector_xcorr_nonref = 0;
  double sum_vector_xcorr = 0;

  //sum the xcorr values of each channel with all others (except itself)
  for(int channel_count1=0; channel_count1< m_numCh; channel_count1++)
  {
    vector_xcorr[channel_count1] = 0;
    for(int channel_count2=0; channel_count2< m_numCh; channel_count2++)
    {
      //do not count itself
      if(channel_count1 != channel_count2)
      {
        vector_xcorr[channel_count1] += m_localXcorr[channel_count1][channel_count2];
      }
    }
  }

	//find the sum of the xcorr values
  for(int channel_count=0; channel_count< m_numCh; channel_count++)
  {
    //only take into account when it is >0
    if(vector_xcorr[channel_count] > 0)
    {
      if(channel_count != (*m_config).reference_channel)
      {
        sum_vector_xcorr_nonref += vector_xcorr[channel_count];
      }
      sum_vector_xcorr += vector_xcorr[channel_count];
    }
  }

  //in weird cases the sum could be 0 (for example when there is no correlation between channels), we avoid that
  //we avoid doing a mean shift to maintain the relation between all samples
  if(sum_vector_xcorr == 0)
	{
		sum_vector_xcorr = 1;
		//all correlations are 0, we force the overall weights to homogeneous
		for(int channel_count=0; channel_count< m_numCh; channel_count++)
			vector_xcorr[channel_count] = 1.0/(double)(m_numCh);
	}

  //convert local xcorrs into weights and adapt the global weight
  float output_weight[m_numCh];
  vector<double> local_xcorr_weight(m_numCh);

  fprintf((*m_fileInOut).featwfd,"0 %d ", frame);
  for(int channel_count=0; channel_count< m_numCh; channel_count++)
  {
    local_xcorr_weight[channel_count] = vector_xcorr[channel_count]/sum_vector_xcorr;

    //output the weights vector into a feature file for further use after the delay&sum
    fprintf((*m_fileInOut).featwfd,"%.4f ",local_xcorr_weight[channel_count]);

    //only adapt when it has not been labelled as noise
    if(!((*m_tdoa).get_delayFilters(channel_count, frame) & F_NOISE))
    {
      m_globalWeight[channel_count] = 0.95 * m_globalWeight[channel_count] + 0.05 * local_xcorr_weight[channel_count];
    }


    //we compute the weight to apply to the output
    if((*m_config).DO_ADAPT_WEIGHTS)
    {
      output_weight[channel_count] = m_globalWeight[channel_count];
    }
    else
    {
      output_weight[channel_count] = 1.0/(float)(m_numCh);
    }

    if((*m_config).DO_AVOID_BAD_FRAMES && (vector_xcorr[channel_count]/sum_vector_xcorr_nonref) < (1.0/((*m_config).badFramesRatio * (m_numCh-1))))
    {
      //we see if there is any channel that has this frame wrong, and eliminate it
        printf("Wrong channel number %d with correlation %f\n", channel_count, (vector_xcorr[channel_count]/sum_vector_xcorr_nonref));
      output_weight[channel_count]=0;
    }
  }
  fprintf((*m_fileInOut).featwfd,"\n");

  //after adaptation we renormalize again (now the global xcorr weight)
  float sum_output_weights = 0;
  for(int channel_count=0; channel_count< m_numCh; channel_count++)
  {
    sum_output_weights += output_weight[channel_count];
    //total_weights += m_outWeightSNR[channel_count];
  }

  //we avoid the case sum=0
  if(sum_output_weights == 0) {sum_output_weights = 1;}

  for(int channel_count=0; channel_count< m_numCh; channel_count++)
  {
      //TEST
      //printf("Weights: %f / %f\n", output_weight[channel_count], sum_output_weights);
    m_channelWeight[channel_count] = output_weight[channel_count]/sum_output_weights; //xcorr weights
  }
}

/*!
  for each pair of channels I compute their xcorr value at the best delay

  \param chan_out Vector containing data from all channels at optimum points
  \return nothing
*/

void DelaySum::compute_local_xcorr_values(vector<vector<float> > &chan_out)
{
  //I suppose that the data will be already in the vectors
  //compute all pairs except the equals

	//first compute the local energy for each channel
	compute_local_energy(chan_out);

  for(int channel_count1 = 0; channel_count1< m_numCh-1; channel_count1++)
  {
    for(int channel_count2 = channel_count1+1; channel_count2<m_numCh; channel_count2++)
    {
      //I use the minimum number of frames among the two channels (in case we are at the end of the file)
      int frames_both;
      if(m_chanFramesRead[channel_count1] > m_chanFramesRead[channel_count2])
      {
        frames_both = m_chanFramesRead[channel_count2];
      }
      else
      {
        frames_both = m_chanFramesRead[channel_count1];
      }

      //compute the xcorr
      m_localXcorr[channel_count1][channel_count2] = 0;
      for(int i=0; i< frames_both; i++)
      {
        m_localXcorr[channel_count1][channel_count2] += chan_out[channel_count1][i]*chan_out[channel_count2][i];
      }
      //normalize by the average energies
      m_localXcorr[channel_count1][channel_count2] /= (m_localEnergy[channel_count1]*m_localEnergy[channel_count2]);


      //don't allow xcorr values to be negative: this means that the correlation is too bad.
      if(m_localXcorr[channel_count1][channel_count2] < 0)
      {
        m_localXcorr[channel_count1][channel_count2] = 0;
      }


      //normalize over all frames
      //m_localXcorr[channel_count1][channel_count2] /= frames_both;
      //symmetrize it
      m_localXcorr[channel_count2][channel_count1] = m_localXcorr[channel_count1][channel_count2];
    }
  }
}

/*!
  for each channel I compute the average energy of the processed window

  \param chan_out Vector containing data from all channels at optimum points
  \return nothing
*/

void DelaySum::compute_local_energy(vector<vector<float> > &chan_out)
{
  for(int channel_count = 0; channel_count< m_numCh; channel_count++)
  {
    m_localEnergy[channel_count] = 0;
    for(int i=0; i< m_chanFramesRead[channel_count]; i++)
    {
			m_localEnergy[channel_count] += chan_out[channel_count][i] * chan_out[channel_count][i];
    }

    //in the weird case that the signal is 0 (only possible if manual editing was done)
    // this signal is considered not interesting by assigning energy = 1
    if(m_localEnergy[channel_count] == 0)
    {
      m_localEnergy[channel_count] = 1;
    }
  }

}

/*!
  Does the sum of the channels with the optimum weights and delays for one segment, overlapping according to the analysis window

  \param best_out 1/2-best TDOa values to output
  \return Nothing
*/

void DelaySum::Channels_Sum(int best_out)
{
    char tmp_string[1024];
  m_log.print_this("Memory allocation for all data\n", 1);
  //for the ICSI shows we need to make sure that any possible skew will not interfere to the
  //memory allocation on the last window
  int max_buffer_size = 2 * m_framesOut + (2 * (*m_config).windowFrames) + m_biggestSkew;

  //vector to hold the output values before writting them to the file, one for each channel
  vector<vector<float> > chan_out;
	chan_out.resize(m_numCh);
  for(int i=0; i<m_numCh; i++)
  {
		chan_out[i].resize(max_buffer_size);
  }

  //single channel floating point buffer
  vector<float> outputf(max_buffer_size);

  //memory for the individual channels output
  vector<vector<float> > outputf_chan;
  for(int channel_count=0; channel_count<m_numCh; channel_count++)
  {
    outputf_chan.push_back(vector<float>(max_buffer_size));
  }

  int frames_to_read;

  int totalNumDelays = (*m_tdoa).get_totalNumDelays();
  m_log.print_this("Finished allocating memory\n", 1);
  //sprintf(tmp_string, "Computed %d delays of %d total\n", totalNumDelays, (*m_tdoa).get_totalNumDelays());fflush(stdout);
  //m_log.print_this(tmp_string, 1);

  //we go through all the computed delays, gathering the data and doing the sum
  for(int count=0; count< totalNumDelays; count++)
  {

    /////initialize the output and move the overlapping segment from the previous window
    for(int i=0; i < m_framesOut; i++)
    {
      outputf[i] = outputf[i+m_framesOut];
      outputf[i+m_framesOut] = 0;
    }
    if((*m_config).INDIV_CHANNELS)
    {
      for(int channel_count=0; channel_count< m_numCh; channel_count++)
      {
        for(int i=0; i < m_framesOut; i++)
        {
          outputf_chan[channel_count][i] = outputf_chan[channel_count][i+m_framesOut];
          outputf_chan[channel_count][i+m_framesOut] = 0;
        }
      }
    }

    //////read the delayed data for each channel
    int data_start_out[MAXNUMCH];
    for(int channel_count=0; channel_count< m_numCh; channel_count++)
    {

      //define the start point to read the data
	  int finalDel = (*m_tdoa).get_finalDelays(channel_count, count, best_out);
      long start_frame = long(count*(*m_config).rate*m_sampleRateInMs) + long(finalDel) + (long)(m_UEMGap) + (long)(m_skew[channel_count]);

      //treat cases at the edges
      data_start_out[channel_count] = 0;
      if(start_frame < 0) //account for the beginning of the file
      {
        data_start_out[channel_count] = -start_frame;
        start_frame = 0;
        for(int i=0; i < data_start_out[channel_count]; i++)
        {
          chan_out[channel_count][i] = 0;
        }
      }

      //define the length to read, taking care of the end of the file
      //we consider the reference channel to find the end of the file
      if(count+1 == totalNumDelays)
      {
        //this is the last block
        frames_to_read = (int)(m_frames-start_frame); //all until the end of each file
      }
      else
      {
        frames_to_read = 2*m_framesOut-data_start_out[channel_count];
      }

      //read the data for each channel
      sf_seek((*m_fileInOut).inputfd[channel_count], start_frame, SEEK_SET);
      m_chanFramesRead[channel_count] = sf_readf_float((*m_fileInOut).inputfd[channel_count], &chan_out[channel_count][data_start_out[channel_count]], frames_to_read);
      if(m_chanFramesRead[channel_count]> max_buffer_size)
      {
        char tmp_string[1024]; sprintf(tmp_string, "Error: Exceeded maximum storage buffer size in frame %d, channel %d! maximum: %d, read: %d\n", count, channel_count, max_buffer_size, m_chanFramesRead[channel_count]);
        m_log.print_this(tmp_string, 10);
        exit(1);
      }

    }

    //we compute the pairwise xcorr of all channels
    compute_local_xcorr_values(chan_out);

    //we compute the weights to be used
    Compute_Sum_Weights(count);


    /////do the weighted sum
    for(int channel_count=0; channel_count< m_numCh; channel_count++)
    {
      float gain_step = m_overallWeight * m_channelWeight[channel_count]/m_framesOut;
      float gain = 0;

      //TEST
      //printf("Weights: %f * %f / %d\n", m_overallWeight, m_channelWeight[channel_count], m_framesOut);

      for(int i=0; i< (m_chanFramesRead[channel_count]+data_start_out[channel_count]); i++)
      {
        outputf[i] += chan_out[channel_count][i] * gain;

        //TEST
        //printf("Output: %f += %f * %f\n", outputf[i], chan_out[channel_count][i], gain);

        //case for the individual channels
        if((*m_config).INDIV_CHANNELS)
        {
          outputf_chan[channel_count][i] += chan_out[channel_count][i] * gain;
        }

        //calculate next gain
        if(i < m_framesOut)
        {
          //going up
          gain += gain_step;
        }
        else if( i >= (m_chanFramesRead[channel_count]+data_start_out[channel_count])-m_framesOut)
        {
          //going down
          gain -= gain_step;
        }
      }
    }

    ////write out the result
    int write_samples_out;
    if(count+1 == totalNumDelays)
    {
      write_samples_out = m_chanFramesRead[(*m_config).reference_channel];
    }
    else
    {
      write_samples_out = m_framesOut;
    }

    //write into the file, will write in whatever format the file is in
    sf_writef_float((*m_fileInOut).outfd[best_out], &outputf[0], write_samples_out);

    //case for the individual channels
    if((*m_config).INDIV_CHANNELS)
    {
      for(int channel_count=0; channel_count< m_numCh; channel_count++)
      {
        sf_writef_float((*m_fileInOut).outputfd[channel_count], &(outputf_chan[channel_count][0]), write_samples_out);
      }
    }

  }

  m_log.print_this("finished writting files\n", 1);

  //force to sync all the outputs
  if((*m_config).INDIV_CHANNELS)
  {
    for(int count=0; count < m_numCh; count++)
    {
      sf_write_sync((*m_fileInOut).outputfd[count]);
    }
  }
  sf_write_sync((*m_fileInOut).outfd[best_out]);

  sprintf(tmp_string, "Sinc-ed all files' outputs for the %d best\n", best_out);
  m_log.print_this(tmp_string, 1);

}

/*!
  Initial channels weights setting algorithm
  \return Nothing

*/

void DelaySum::Set_Channels_Weights()
{
    //I initialize it
    for(int channel_count=0; channel_count< m_numCh; channel_count++)
        m_globalWeight[channel_count] = 1/(float)(m_numCh);
}

/*
  checks that we don't overpass the end of the file when accessing a particular data frame,
  \param counter Segment number to be processed
  \return the start point or -1 if before start

*/

long DelaySum::get_start_frame(long counter)
{
  //checks that we don't overpass the end of the file, and returns the start point or -1
  long start_frame = (long)(counter*(*m_config).rate*m_sampleRateInMs + m_UEMGap);

  if((start_frame + (*m_config).windowFrames + m_biggestSkew) < m_frames)
  {
    //we are ok
    return start_frame;
  }
  else
  {
    return -1;
  }

}

/*!
  Wrapper for the get_channels_data when no delays are given, therefore delays 0
  \param file_start frame where to start retrieval of the channels data
  \return nothing
*/
void DelaySum::get_channels_data(long file_start)
{
  vector<int> tmp_delays(m_numCh, 0);
  get_channels_data(file_start, tmp_delays);
}


/*!
  Retrieve the acoustic data from the acoustic files for all channels
  \param file_start frame where to start retrieval of the channels data
  \param ch_delays delay to be applied to each channel when retrieving the data
  \return Nothing
*/

void DelaySum::get_channels_data(long file_start, vector<int> ch_delays)
{
  long vector_start = 0; //position in the vector where to start writting data from the file
  int real_window; // number of frames that are read from the file. It is = window_size if there is no skew
  int skew;

  //file start is where to start to read from the file

  for(int count=0; count< m_numCh; count++)
  {
    //real window is different from window_size in case we have a negative skew
    real_window = (*m_config).windowFrames;
    vector_start = 0;
    skew = m_skew[count];
    //we consider the possible skew of the channels in the ICSI meeting. For other meetings the skew is 0
    if((file_start + ch_delays[count] + skew)<0)
    {
      //some frames at the beginning fall before the beginning of the file.
      vector_start = -(file_start + ch_delays[count] + skew);
      file_start = - (skew + ch_delays[count]);
      real_window += (skew + ch_delays[count]);
    }

    //if((*m_fileInOut).readChanData(count, &m_chanData[count][vector_start], file_start+skew+ch_delays[count], real_window) != real_window)
    int readCount = (*m_fileInOut).readChanData(count, &m_chanData[count][vector_start], file_start+skew+ch_delays[count], real_window);
    if( readCount != real_window)
    {
      char tmp_string[1024]; sprintf(tmp_string,"Error retrieving data for channel %d, not enough samples in the file\n",count);
      m_log.print_this(tmp_string, 10);
      exit(1);
    }

    //assign it to the correct position
    if(vector_start != 0)
    {
      //we fill the initial elements with 0'os
      for(int i=0; i<vector_start; i++)
      {
        m_chanData[count][i] = 0;
      }
    }
    //save the number of frames read
    m_chanFramesRead[count] = real_window;
  }
}

/*!
  Starting function to compute the TDOA values for a given position
  \param start_frame starting frame to compute the TDOA values
  \param delayNum number of Nbest delays to compute
  \return nothing
*/
void DelaySum::computeTDOAValues(long start_frame, int delayNum)
{
    ///we retrieve the data for all the channels
    get_channels_data(start_frame);

    ///we do the xcorrelation, it returns N maximum delays in order from max to min. We don't do it for the reference one, it's set to 0

    vector<int> bestDelays = (*m_tdoa).xcorrelation_FFTReal_full(delayNum, m_chanData);

    //recompute the reference channel given the selected delays
    if((*m_config).COMPUTE_REFERENCE == 2)
    {
      (*m_config).reference_channel = find_best_channel_xcorr_adapt(start_frame, bestDelays);
    }
}

/*!
	Evaluation function that computes the difference between the two imput channels and writes it to the output file.
	It can be used to compare 2 different runs of beamformit.

*/
void DelaySum::computeChannelsDifference()
{
  int max_buffer_size = 2 * m_framesOut;
	float signalPower = 0;

	if(m_numCh != 2)
	{
		char tmp_string[1024];
		sprintf(tmp_string, "Number of channels (%d) is different than 2, incompatible with this mode\n", m_numCh);
	  m_log.print_this(tmp_string, 10);
    exit(1);
	}

  //vector to hold the output values before writting them to the file, one for each channel
  vector<vector<float> > chan_out;
  for(int i=0; i<m_numCh; i++)
  {
    chan_out.push_back(vector<float>(max_buffer_size));
  }

  //single channel floating point buffer
  vector<float> outputf(max_buffer_size);
  vector<float> outputf2(max_buffer_size);

  int frames_to_read;
  for(int frameCount=0; frameCount < m_frames; frameCount += m_framesOut)
  {
		if(frameCount + m_framesOut < m_frames)
			frames_to_read = m_framesOut;
		else
			frames_to_read = m_frames - frameCount; //last block

    //////read the data for each channel
    for(int channel_count=0; channel_count< m_numCh; channel_count++)
    {
      //read the data for each channel
      sf_seek((*m_fileInOut).inputfd[channel_count], frameCount, SEEK_SET);
      m_chanFramesRead[channel_count] = sf_readf_float((*m_fileInOut).inputfd[channel_count], &chan_out[channel_count][0], frames_to_read);
      if(m_chanFramesRead[channel_count]> max_buffer_size)
      {
        char tmp_string[1024];
				sprintf(tmp_string, "Error: Exceeded maximum storage buffer size in frame %d, channel %d! maximum: %d, read: %d\n", frameCount, channel_count, max_buffer_size, m_chanFramesRead[channel_count]);
        m_log.print_this(tmp_string, 10);
        exit(1);
      }
    }

		//compute the output signal
		int length;
		if(m_chanFramesRead[0] > m_chanFramesRead[1])
			length = m_chanFramesRead[1];
		else
			length = m_chanFramesRead[0];
		for(int i=0; i < length; i++)
		{
			outputf[i] = chan_out[0][i] - chan_out[1][i];
			outputf2[i] = chan_out[0][i] / chan_out[1][i];
			signalPower += outputf[i] * outputf[i];
		}

    ////write out the result
    //write into the file, will write in whatever format the file is in
	sf_writef_float((*m_fileInOut).outfd[0], &outputf[0], length);
	sf_writef_float((*m_fileInOut).outfd[1], &outputf2[0], length);
  }

  //force to sync
  sf_write_sync((*m_fileInOut).outfd[0]);
  sf_write_sync((*m_fileInOut).outfd[1]);

	//output the total power of the residual
  char tmp_string[1024]; sprintf(tmp_string, "Residual average power: %f\n", signalPower/m_frames);
  m_log.print_this(tmp_string, 10);

}
