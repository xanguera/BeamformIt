/*
	Holds the TDOA and XCORR data and the set of filters to apply after the TDOA computation 
	in order to stabilize them and select the output delay
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>

#include "tdoa.h"
#include "fileinout.h"
#include "global.h"

/*!
  Constructor
  \param fileInOut in/out class
  \param config configuration class
*/
TDOA::TDOA(FileInOut* fileInOut, Configuration* config)
{
	m_fileInOut = fileInOut;
	m_config = config;
}

/*!
	Initializes the structures to be used with the class

	\param UEMGap UEM gap
	\return Nothing
*/
void TDOA::init(int numChannels, float UEMGap)
{
	//initialize the passed values
	m_numCh = numChannels;
	m_UEMGap = UEMGap;

	//create dummy vectors
	vector<int> dummy_int_nbest((*m_config).nbest_amount,0);
	vector<float> dummy_float_nbest((*m_config).nbest_amount,0);
	dummy_float_nbest[0]=1; //by default the first value is the best
	vector<int> dummy_int_2best(2,0);
	vector<float> dummy_float_2best(2,0);
	vector<float> dummy_float_dels(m_totalNumDelays,0);

	vector<unsigned int> dummy_Uint_dels(m_totalNumDelays,0);
	vector<vector<int> > dummy_int2(m_totalNumDelays, dummy_int_nbest);
	vector<vector<float> > dummy_float2(m_totalNumDelays, dummy_float_nbest);
	vector<vector<int> > dummy_int3(m_totalNumDelays, dummy_int_2best); //for the 2-best
	vector<vector<float> > dummy_float3(m_totalNumDelays, dummy_float_2best);
	vector<int> dummy_int_dels(m_totalNumDelays,0);


	m_bestClass.resize(m_numCh, dummy_Uint_dels);
	m_bestClass2.resize(m_numCh, dummy_Uint_dels);

	//then for the channel level
  
	m_chDelays.resize(m_numCh, dummy_int2);
	m_chXcorrValues.resize(m_numCh, dummy_float2);
	m_delayFilters.resize(m_numCh, dummy_int_dels);
	m_finalDelays.resize(m_numCh, dummy_int3);
	m_finalXcorrValues.resize(m_numCh, dummy_float3);
        m_globalDelayFilters.resize(m_totalNumDelays, 0);

  //I find a window for the fft that is the inmediate upper power of 2 from the window parameter
  m_fftWindowFrames = int(pow(2.0, ceil(log(double((*m_config).windowFrames))/log(2.0))));
  char tmp_string[1024]; sprintf(tmp_string, "Using an fft window of %d points for the window of %d points\n", m_fftWindowFrames, (*m_config).windowFrames);
  m_log.print_this(tmp_string, 1); 

}

/*!
  we do all the possible filterings to enhance the delays
  each of the filters alters the computed delays and xcorr_values
  Only the vector m_bestClass is modified with the 
  selected delay, except for the noise_filter, which modifies the source vector.
  
  \return Nothing
*/

void TDOA::Optimum_delays_filtering()
{

  float threshold;
  
  //We consider two kinds of noise thresholding, a variable threshold at 10% of the
  // possible xcor values, a fixed one to 0.1 or none at all.

  //level out the noise parts by assigning thse samples to the previous element
  if((*m_config).DO_NOISE_THRESHOLD == 1)
  {
    threshold = determine_noise_threshold();
    char tmp_string[1024]; sprintf(tmp_string, "Thresholding noisy frames lower than %f\n", threshold);
    m_log.print_this(tmp_string, 1);
  }
  else if((*m_config).DO_NOISE_THRESHOLD == 2)
  {
    threshold = (*m_config).min_xcorr;
  }
  else
  {
    m_log.print_this("Not applying any noise thresholding\n", 1);
    return;
  }

  m_log.print_this("  Aplying a noise thresholding\n", 1);
  noise_filter(threshold);
  m_log.print_this("  Finished applying a noise filter\n", 1);
}


/*!
  Determines the threshold to apply to avoid noisy frames
  \return Threshold value
*/

float TDOA::determine_noise_threshold()
{
  // we define a global threshold by averaging all dimensions of the delays vector
  //first define the thresholds for each channel
  float threshold;
  vector<float> sorted_values(m_totalNumDelays,0);
  for(int channel_count=0; channel_count< m_numCh; channel_count++)
  {
    //assign all 1st-best values
    for(int frame=0; frame<m_totalNumDelays; frame++)
      if(m_chXcorrValues[channel_count][frame][m_bestClass[channel_count][frame]] != 1)
        sorted_values[frame] += m_chXcorrValues[channel_count][frame][m_bestClass[channel_count][frame]];
  }  
  //sort the vector and select the threshold
  sort(sorted_values.begin(), sorted_values.end());
  threshold = sorted_values[(int)((*m_config).noise_percent*m_totalNumDelays/100)]/(m_numCh-1);
  char tmp_string[1024]; sprintf(tmp_string, "Threshold is %f (min: %f max %f)\n", threshold, sorted_values[0]/(m_numCh-1), sorted_values[m_totalNumDelays-1]/(m_numCh-1));
  m_log.print_this(tmp_string, 1);

  return(threshold);

}

/*!
  Applies a noise filter based on thresholding the xcorr values to avoid using non reliable TDOa data
  \param threshold Threshold to be applied, from 0 to 1
*/

void TDOA::noise_filter(float threshold)
{
  //apply a minimum value for the xcorrelation. Assign the previous delay in case it doesn't pass

  for(int channel_count=0; channel_count < m_numCh; channel_count++)
  {
    //go through all computed delays
    for(int count=0; count<(int)(m_chDelays[channel_count].size()); count++)
    {
	  //check wether the xcorr value for the best match is below the threshold
	  // when it is, I propagate the previous values forward, smoothing the delays
      if(m_chXcorrValues[channel_count][count][0] < threshold)
      {
		//printf("Channel %d Frame %d, Xcorr %f is below threshold %f\n", channel_count, count, m_chXcorrValues[channel_count][count][0], threshold);
        //enable the flag saying that a noise threshold hs been applied
        m_delayFilters[channel_count][count] += F_NOISE;
	      
        // I don't repeat the value as I let the viterbi handle the situation
        if(count>0)
        {
          //I just copy the delay value, but not the xcorr value
          m_chDelays[channel_count][count] = m_chDelays[channel_count][count-1];
        }
        else
        {
          //if the first delay is silence, we turn it to 0;
          m_chDelays[channel_count][count].assign(m_chDelays[channel_count][count].size(), (*m_config).marginFrames); //a big discouraging value (the max we allow as delay)
          m_chDelays[channel_count][count][0]=0;
          m_chXcorrValues[channel_count][count].assign(m_chXcorrValues[channel_count][count].size(),0);
          m_chXcorrValues[channel_count][count][0]=1;
        }     
      }
	//printf("Channel %d frame %d: %d %d %d %d\n", channel_count, count, m_chDelays[channel_count][count][0], m_chDelays[channel_count][count][1], m_chDelays[channel_count][count][2], m_chDelays[channel_count][count][3]);
    }
  }
}

/*******************************************************************************
* Double Viterbi acoustic Modeling
********************************************************************************/

/*!
  does the full delays modelling using the cluster program's routines
  
  **NOT ANYMORE, given variable reference channel**
  we eliminate the reference channel's delays, which are always the same and don't need modelling
  NEED to be consistent with this in the Viterbi

  \return Nothing
*/

void TDOA::Optimum_delays_acoustic_modelling()
{

  //we start with an Nbest individual channels viterbi, which will detect if there is any
  // delay that got pushed back to 3rd or more order during detection.
  //as this is channel independent, it doesn't require the computation that the other one does
  Nbest_IndivCh_decoding();

  //In this second part we do a 2-best (better due to computation requirements) of the delays on all channels
  Global_delays_decoding();

  //perform an overlap detection: EXPERIMENTAL CODE
  if((*m_config).DO_OUTPUT_OVERLAP)
      Compute_overlap_models();
  
  m_log.print_this("Finished computing optimum delays\n", 1);

}

/*!
  does an Nbest decoding of each individual channel in order to reduce the amount of 
  computation of the full decoding without loosing information from all channels.

  \return nothing
  
*/

void TDOA::Nbest_IndivCh_decoding()
{
  m_log.print_this("Deconding Each channel Independently to find the 2-best delays from each N-best input", 1);fflush(stdout);
  segment_Nbestdelays_2best();
  m_log.print_this("Finished decoding both best paths, writting out results\n", 1);fflush(stdout);

  //now we need to modify the delays vector to put the selected delays within the 2-best for the next processing
  vector<int> best_del((*m_config).nbest_amount,0);
  vector<float> best_val((*m_config).nbest_amount,0);
  unsigned int position;
  for(int frame=0; frame<m_totalNumDelays; frame++)
  {
    for(int channel_count=0; channel_count< m_numCh; channel_count++)
    {
      //currently need to compute for all, as ref. channel can vary
      //if(channel_count != (*m_config).reference_channel)
      //{
        //we copy at the beginning the 2-best and then the others following
        best_del[0] = m_chDelays[channel_count][frame][m_bestClass[channel_count][frame]];
        best_del[1] = m_chDelays[channel_count][frame][m_bestClass2[channel_count][frame]];
        best_val[0] = m_chXcorrValues[channel_count][frame][m_bestClass[channel_count][frame]];
        best_val[1] = m_chXcorrValues[channel_count][frame][m_bestClass2[channel_count][frame]];
        position = 2;	      
        for(unsigned int i=0; i<(*m_config).nbest_amount; i++)
        {
          if(i != m_bestClass[channel_count][frame] && i != m_bestClass2[channel_count][frame] && position < (*m_config).nbest_amount)
          {
            best_del[position] = m_chDelays[channel_count][frame][i];
            best_val[position] = m_chXcorrValues[channel_count][frame][i];
            position++;
          }
        }
	      
        //we copy everything back to place
        for(unsigned int i=0; i<(*m_config).nbest_amount; i++)
        {
          m_chDelays[channel_count][frame][i] = best_del[i];
          m_chXcorrValues[channel_count][frame][i] = best_val[i];
        }
      //}
      //and we also reset the best delay
      m_bestClass[channel_count][frame] = 0;
    }
  }
}

/*!
  Second pass multichannel Viterbi decoding of the 2-best individual channels TDOA values combinations

  \return nothing

*/

void TDOA::Global_delays_decoding()
{
  m_log.print_this("\nDecoding all 2-best input delays to find the optimum path in coherence with all channels\n", 1);

  //maximum size of the delays vector
  int max_block_size = 5;

  //There are some conditions where there is a HUGE amount of microphones available for processing.
  // In such cases we need to do a suboptimal processing of the 2-best cases by doing it by blocks
  //where each block uses the already optimized values from the previous block.
  int num_channels = m_numCh; //changed from m_numCh-1 as we now consider also the reference
  if(num_channels > max_block_size)
  {
    //first we search for the optimum chunk size (to minimize computation)
    int optimum_block_size = (int)((float)(num_channels)/((int)(floor((float)(num_channels)/max_block_size)) + 1));
    char tmp_string[1024]; sprintf(tmp_string, "Optimum block size selected at %d\n", optimum_block_size);
    m_log.print_this(tmp_string, 1);
    
    int offset = 0, block_size;
    while(offset <(num_channels))
    {
      if(offset+2*optimum_block_size > (num_channels))
      {
        block_size = (num_channels) - offset;
      }
      else
      {
        block_size = optimum_block_size;
      }
      char tmp_string[1024]; sprintf(tmp_string, "Block-based Viterbi segmentation using channels %d to %d\n", offset, offset+block_size);
      m_log.print_this(tmp_string, 1);
      segment_Nbestdelays_multichannel(2, block_size, offset);
      offset += block_size;
    }	  
  }
  else
  {
    //all in one block
    segment_Nbestdelays_multichannel(2, num_channels, 0);
  }  
}


/*! 
  Performs an overlap detection based on several criteria, outputs results into a log file

  \return Nothing

*/

void TDOA::Compute_overlap_models()
{
  int window = (int)(500.0/(*m_config).rate); //half window in ms
  int nSamples = m_totalNumDelays;
  vector<float> xcorr_value1(nSamples, 0); // aritmetic mean
  vector<float> xcorr_value2(nSamples, 0); // geometric mean
  vector<float> xcorr_value3(nSamples, 0); // entropy

  //we compute the values for all frames
  for(int frame=0; frame < nSamples; frame++)
  {
    xcorr_value1[frame] = 0;
    xcorr_value2[frame] = 0;

    for(int channel_count=0; channel_count< m_numCh; channel_count++)
    {
      if(channel_count != (*m_config).reference_channel)
      {
        if(m_bestClass[channel_count][frame] == 0)
        {
          xcorr_value1[frame] += m_chXcorrValues[channel_count][frame][1];
          xcorr_value2[frame] += log(m_chXcorrValues[channel_count][frame][1]);
          if(m_chXcorrValues[channel_count][frame][1] != 0)
          {
            xcorr_value3[frame] += m_chXcorrValues[channel_count][frame][1]*log(m_chXcorrValues[channel_count][frame][1]);
          }
        }
        else
        {
          xcorr_value1[frame] += m_chXcorrValues[channel_count][frame][0];
          xcorr_value2[frame] += log(m_chXcorrValues[channel_count][frame][0]);		  
          if(m_chXcorrValues[channel_count][frame][0] != 0)
          {
            xcorr_value3[frame] += m_chXcorrValues[channel_count][frame][0]*log(m_chXcorrValues[channel_count][frame][0]);
          }
        }
      }	  
    } 
    xcorr_value1[frame] /= (m_numCh-1);
    xcorr_value2[frame] = exp(xcorr_value2[frame]/(m_numCh-1));
    xcorr_value3[frame] /= -(m_numCh-1);
  }
  
  vector<float> temp_xcorr = xcorr_value1;
  sort(temp_xcorr.begin(), temp_xcorr.end());
  float threshold = temp_xcorr[(int)(temp_xcorr.size()*0.90)];
  char tmp_string[1024]; sprintf(tmp_string, "Xcorr at 90 percent : %f\n", temp_xcorr[(int)(temp_xcorr.size()*0.90)]);
  m_log.print_this(tmp_string, 1);


  //compute the average over a window and define if it is overlap or not
  vector<int> overlap(nSamples, 0);
  //float threshold = 0.1;
  float ave_xcorr1;
  float ave_xcorr2;
  float ave_xcorr3;
  for(int frame=0; frame < nSamples; frame++)
  {
    int start = frame - window; if(start<0){start=0;}
    int end = frame + window; if(end>nSamples){end = nSamples;}
    ave_xcorr1 = ave_xcorr2 = ave_xcorr3 = 0;
    for(int i=start; i<end; i++)
    {
      ave_xcorr1 += xcorr_value1[i];
      ave_xcorr2 += xcorr_value2[i];
      ave_xcorr3 += xcorr_value3[i];
    }
    ave_xcorr1 /= (end-start);
    ave_xcorr2 /= (end-start);
    ave_xcorr3 /= (end-start);
    if(ave_xcorr1 > threshold)
    {
      overlap[frame] = 1;
    }
    char tmp_string[1024]; sprintf(tmp_string, " %d", overlap[frame]);
    m_log.print_this(tmp_string, 0);
  }

  //expand the overlap vector to represent the status of each frame in absolute value (independent of the
  // window used): by doing so we filter also spurious peaks
  
  vector<int> overlap_out(nSamples, 0);
  for(int frame = 0; frame < nSamples; frame++)
  {
    if(overlap[frame] == 1)
    {
      //avoid spurious
      if((overlap[frame-1] == 1 && frame > 0) || (overlap[frame+1] == 1 && frame<nSamples-1))
      {
        int start = frame - window; if(start<0){start = 0;}
        int end = frame + window; if(end>nSamples){end = nSamples;}
        for(int j=start; j<end; j++)
        {
          overlap_out[j] = 1;
        }
      }
    }
  }  
  m_log.print_this("============\n", 0);

  //output an rttm-like file with all overlap regions
  int start_ovl = 0;
  int end_ovl = 0;
  for(int frame=1; frame < nSamples; frame++)
  {
    char tmp_string[1024]; sprintf(tmp_string, " %d", overlap_out[frame]);
    m_log.print_this(tmp_string, 0);
    //0 to 1
    if(overlap_out[frame-1] == 0 && overlap_out[frame] == 1)
    {
      start_ovl = frame;
    }
    //1 to 0
    if(overlap_out[frame-1] == 1 && overlap_out[frame] == 0)
    {
      end_ovl = frame;
      fprintf((*m_fileInOut).ovlfd, "%s 1 %.2f %.2f\n", (*m_config).SHOWNAME.c_str(), ((float)((*m_config).rate*start_ovl)/1000), (float)(end_ovl*(*m_config).rate)/1000);
    }
  }  
  //the last one
  if(overlap_out[nSamples-1] == 1)
  {
    end_ovl = nSamples-1;
    fprintf((*m_fileInOut).ovlfd, "%s 1 %.2f %.2f\n", (*m_config).SHOWNAME.c_str(), ((float)((*m_config).rate*start_ovl)/1000), (float)(end_ovl*(*m_config).rate)/1000);
  }

}

/*!
  This is a Viterbi decoding of the Nbest delays for each individual channel
  the probability of each delay is taken from the GCC-PHAT xcorr value, normalized to
  summ up 1. The transition probability is a function of the distance between 2 
  adjacent delays.
  we have config.nbest_amount models
  
  \return Nothing

*/

void TDOA::segment_Nbestdelays()
{

  //data structures declaration and memory allocation
  int nSamples = m_totalNumDelays;
  int nClass = (*m_config).nbest_amount;
  int min_duration = 1;
  float trans_weight = (*m_config).trans_weight_nbest; //weight to apply to the transition probability
  char tmp_string[1024]; sprintf(tmp_string, "\nTrans_weight selected to %f\n", trans_weight);
  m_log.print_this(tmp_string, 1);
  int totStates = min_duration*nClass; 
  vector<int> F(nSamples,0); // hold the frame backpointer
  vector<int> R(nSamples,0); // hold the best cluster at each frame
  vector<vector<double> > dC;// hold the best scores for each state
  dC.push_back(vector<double>(totStates, 0));
  dC.push_back(vector<double>(totStates, 0));
  vector<vector<int> > tC;// hold the frame pointers for each state
  tC.push_back(vector<int>(totStates, 0));
  tC.push_back(vector<int>(totStates, 0));
  vector<double> lkld(nClass,0);  //lkld probabilities
  vector<vector<double> > trans; //transition probabilities from all previous to all current
  for(int i=0; i<nClass; i++)
  {
    trans.push_back(vector<double>(nClass,0));
  }

  sprintf(tmp_string, "Starting Viterbi with %d channels and %d N-best values per channel\n", m_numCh, nClass);
  m_log.print_this(tmp_string, 1);

  for(int channel_count=0; channel_count < m_numCh; channel_count++)
  {
    //if(channel_count != (*m_config).reference_channel)
    //{
      char tmp_string[1024]; sprintf(tmp_string, "Computing Viterbi for channel %d\n", channel_count);
      m_log.print_this(tmp_string, 1);
	      
      F[0] = 0;
      R[0] = 0;
	  
      // first initialize the global scores for all states to a small likelihood
      for(int st=0; st<totStates; st++){
        dC[LAST][st] = -1000;
        dC[CURR][st] = -1000;
        tC[LAST][st] = 0;
        tC[CURR][st] = 0;
      }
	  
      //now set the lkld for the 0th frame and set the 0th state
      //we normalize the xcorr values
      //in here there is no transition prob, it is the start.
      lkld = compute_xcorr_prob(0, channel_count);
	  	  
      for(int cl=0; cl<nClass;cl++)
      {
        dC[LAST][(cl*min_duration)] = lkld[cl];
        dC[CURR][(cl*min_duration)] = lkld[cl];
      }
	  
      // 
      // Loop over remaining frames of data
      //
      for(int frame = 1; frame <nSamples; frame++) {

        //compute lkld.
        lkld = compute_xcorr_prob(frame, channel_count);	    

        //we need the trans. prob from all states in the previous frame to al states
	//printf("Channel %d\n", channel_count);        
	trans = compute_trans_prob(frame, channel_count);

        // the index of the final state of the best previous cluster
        int bestGlobalState = R[frame-1]*min_duration + min_duration-1;

        // 
        // Handle state 0 of each model
        //   This state *must* come from the final state of best model
        //   at the previous time.
        // It only makes sense if we have a minimum duration>1
        //
        if(min_duration > 1)
        {
          for(int cl=0; cl<nClass; cl++){
            int curState;
            // we compute the transition prob for that cluster from the previous best
            // to the current computed.
            double clusterTrans = dC[LAST][bestGlobalState] + trans_weight * trans[R[frame-1]][cl];

            //
            // Get the index of the 0th state of current cluster
            //
            curState = cl*min_duration;	
		  
            // add local distance to state's global score
            dC[CURR][curState] = clusterTrans + (double)(lkld[cl]);
            // since we are transitioning now, this frame is the begin time
            tC[CURR][curState] = frame;
          }
        }

        // 
        // Handle the other states of each cluster
        //
        for(int cl=0; cl<nClass;cl++) {
          double forwardTrans;
          int gl_state;		// index into global state array
          int cl_state;		// index of state in a cluster
          double localDist = (double)(lkld[cl]);
	      
          //
          // Start with state 1 and end at the second-to-last state
          //   There are no self-loops on the internal states!
          //
          for(cl_state = 1; cl_state < (min_duration-1); cl_state++) 
          {
            gl_state = (cl*min_duration) + cl_state;
            //WARNING: We consider a trans probability so that it counts in each
            // step of the process, regardless of the min. duration.
            //it could be understood as an addendum to the state lkld
            forwardTrans = dC[LAST][gl_state-1] + trans_weight * trans[cl][cl]; 
            dC[CURR][gl_state] = forwardTrans + localDist;
            tC[CURR][gl_state] = tC[LAST][gl_state-1];
          }
	      
          //
          // Now do the final state of each cluster
          //  This is handled differently from the internal states because
          //  it has a self-loop prob
          // this changes slightly if min. dur. is 1
          //
          if(min_duration > 1)
          {
            gl_state = (cl*min_duration) + cl_state;
            forwardTrans = dC[LAST][gl_state-1] + trans_weight * trans[cl][cl];
          }
          else
          {
            //I come from the best in the prev run
            gl_state = cl;
            forwardTrans = dC[LAST][bestGlobalState] + trans_weight * trans[R[frame-1]][cl];
          }
          double selfLoopTrans = dC[LAST][gl_state] + trans_weight * trans[cl][cl];

          if (selfLoopTrans >= forwardTrans) {
            dC[CURR][gl_state] = selfLoopTrans + localDist;
            tC[CURR][gl_state] = tC[LAST][gl_state];
          } else {
            dC[CURR][gl_state] = forwardTrans + localDist;
            if(min_duration>1)
            {
              tC[CURR][gl_state] = tC[LAST][gl_state-1];
            }
            else
            {
              //we transition now, we point to here as the start of the cluster
              tC[CURR][gl_state] = frame;
            }
          }
        }

        // 
        // Finished with all states of all models. Now find the cluster
        // that has the highest likelihood in its final state.
        //
        double MaxLkld = -1e30;
        int MaxClass = -1;
        int MaxClassFinSt = -1;
        for(int cl=0; cl<nClass; cl++) {
          int finalClusState = cl*min_duration + min_duration-1;
          double tmplkld = dC[CURR][finalClusState];
          if(tmplkld > MaxLkld) {
            MaxClass = cl;
            MaxLkld = tmplkld;
            MaxClassFinSt = finalClusState;
          }
        }
	    
        if (MaxClass == -1) {
          m_log.print_this("ERROR: error finding best class\n", 10);
          exit(1);
        }
        // Record the results
        R[frame] = MaxClass;
        F[frame] = tC[CURR][MaxClassFinSt];

        // copy CURR to LAST
        for(int st=0;st <totStates; st++) {
          dC[LAST][st] = dC[CURR][st];
          tC[LAST][st] = tC[CURR][st];
        }

      } /* frame loop */


      //last one
      int from = F[nSamples-1];
      int to;
      for(int i = from; i<nSamples; i++)
      {
        m_bestClass[channel_count][i] = R[nSamples-1];
	//m_bestClass[channel_count][i] = m_chDelays[channel_count][i][R[nSamples - 1]];
      }

      //others
      while(from>1)
      {
        to = from-1; //included
        from = F[to];
        for(int i=from; i<=to; i++)
        {
          m_bestClass[channel_count][i] = R[to];
          //m_bestClass[channel_count][i] = m_chDelays[channel_count][i][R[to]];
        }
      }

      //print out some results
      int bestFinalClus = R[nSamples-1];
      int bestFinalState = bestFinalClus * min_duration + min_duration - 1;
      double viterbi_score = (dC[CURR][bestFinalState]); //I cast back from double to float
      sprintf(tmp_string, "   Segmentation Over with Viterbi score = %f\n", viterbi_score);
      m_log.print_this(tmp_string, 1);
      fflush(stdout);
	
    //}
  }
}

/*!
  Computes the Viterbi in the individual channels twice, to obtain the 2-best delays for each channel independently

  \return Nothing

*/

void TDOA::segment_Nbestdelays_2best()
{
  
  //first we run the normal way to get the 1st best
  m_log.print_this("Computing the 1st best delay\n", 1);
  printf("First best start\n");fflush(stdout);
  segment_Nbestdelays();
  printf("First best finish\n");fflush(stdout);

  //now we copy the chosen delays into a tmp vector and their xcorr values
  std::vector<std::vector<unsigned int> > best_class_tmp;
  vector<unsigned int> dummy_Uint_dels(m_totalNumDelays,0);
  best_class_tmp.resize(m_numCh, dummy_Uint_dels);

  std::vector<std::vector<float> > best_xcorr_tmp;
  vector<float> dummy_float_dels(m_totalNumDelays,0);
  best_xcorr_tmp.resize(m_numCh, dummy_float_dels);

  //now we discourage the selected values from being decoded again
  for(int channel_count = 0; channel_count<m_numCh; channel_count++)
  {
    for(int frame=0; frame<m_totalNumDelays; frame++)
    {
      best_xcorr_tmp[channel_count][frame] = m_chXcorrValues[channel_count][frame][m_bestClass[channel_count][frame]];
      best_class_tmp[channel_count][frame] = m_bestClass[channel_count][frame];
      //test
      //printf("1step Viterbi Selected channel %d frame %d: %d (%f)\n", channel_count, frame, m_chDelays[channel_count][frame][best_class_tmp[channel_count][frame]], best_xcorr_tmp[channel_count][frame]);fflush(stdout);
      m_chXcorrValues[channel_count][frame][m_bestClass[channel_count][frame]] = 0;
    }
  }
  
  //now we run the decoding again, hoping for a new set of best values
  m_log.print_this("Computing the 2nd best delay\n", 1);fflush(stdout);
  printf("second best start\n");fflush(stdout);
  segment_Nbestdelays();
  printf("second best finish\n");fflush(stdout);
  m_log.print_this("Done computing the 2nd best delay\n", 1);fflush(stdout);

  //finally, we copy all the values to the appropriate places and restore the 1st-best xcorr values
  for(int channel_count = 0; channel_count<m_numCh; channel_count++)
  {
    for(int frame=0; frame<m_totalNumDelays; frame++)
    {
      m_bestClass2[channel_count][frame] = m_bestClass[channel_count][frame];
      m_bestClass[channel_count][frame] = best_class_tmp[channel_count][frame];
      m_chXcorrValues[channel_count][frame][m_bestClass[channel_count][frame]] = best_xcorr_tmp[channel_count][frame];
    }
  }  
}

/*!
  Runs a viterbi segmentation taking into account all channels at the same time and 
  using a gaussian model to model the delays. We pass as a parameter what is the Nbest to treat
  and 2 parameters defining a block size and an offset to reduce the computation when I have a lot of channels

  This is a Viterbi decoding of the Nbest delays for each individual channel.
  The transition probability is a function of the distance between 2 
  adjacent delays.
  We consider all possible combinations between the N-best values of all channels except the 
  reference channel

  \param Nbest Number of best channels to compute
  \param block_size number of channels involved in the Viterbi
  \param offset offset from channel 0 to begin counting the channels to process
  \return Nothing

*/

void TDOA::segment_Nbestdelays_multichannel(int Nbest, int block_size, int offset)
{

  //data structures declaration and memory allocation
  int nSamples = m_totalNumDelays;
  //the possible classes are now all possible combinations of the Nbest values from all mics except the reference one
  //int nClass = (int)(pow((double)(Nbest), (double)(m_numCh-1)));
  int nClass = (int)(pow((double)(Nbest), (double)(block_size)));
  char tmp_string[1024]; sprintf(tmp_string, "Total number of possible delays: %d\n", nClass);
  m_log.print_this(tmp_string, 1);

  int min_duration = 1;
  float trans_weight = (*m_config).trans_weight_multi; //weight to apply to the transition probability
  int totStates = min_duration*nClass; 
  vector<int> F(nSamples,0); // hold the frame backpointer
  vector<int> R(nSamples,0); // hold the best cluster at each frame
  vector<vector<double> > dC;// hold the best scores for each state
  dC.push_back(vector<double>(totStates, 0));
  dC.push_back(vector<double>(totStates, 0));
  vector<vector<int> > tC;// hold the frame pointers for each state
  tC.push_back(vector<int>(totStates, 0));
  tC.push_back(vector<int>(totStates, 0));
  //lkld probabilities for all possible combinations
  vector<double> lkld(nClass,0);  //lkld probabilities
  //transition probabilities from all previous to all current for ALL channels (in vector format)
  //first element is "from" and second is "to"
  vector<vector<double> > trans; 
  for(int i=0; i<nClass; i++)
  {
    trans.push_back(vector<double>(nClass,0));
  }

  //vector with all possible N-best combinations
  vector<vector<int> > combinations;


  m_log.print_this("Starting the viterbi\n",1);

  F[0] = 0;
  R[0] = 0;
  
  // first initialize the global scores for all states to a small likelihood
  for(int st=0; st<totStates; st++){
    dC[LAST][st] = -1000;
    dC[CURR][st] = -1000;
    tC[LAST][st] = 0;
    tC[CURR][st] = 0;
  }
  
  
  //now set the lkld for the 0th frame and set the 0th state
  //we normalize the xcorr values
  //in here there is no transition prob, it is the start.
  combinations = create_delays_positions_block(Nbest, block_size, offset, 0);
  lkld = compute_xcorr_prob_multichannel(0, combinations);
  
  
  for(int cl=0; cl<nClass;cl++)
  {
    dC[LAST][(cl*min_duration)] = lkld[cl];
    dC[CURR][(cl*min_duration)] = lkld[cl];
  }
  
  // 
  // Loop over remaining frames of data
  //
  for(int frame = 1; frame <nSamples; frame++) {
    
    //create the delays positions for that frame
    combinations = create_delays_positions_block(Nbest, block_size, offset, frame);
    //compute lkld.
    lkld = compute_xcorr_prob_multichannel(frame, combinations);

    //we need the trans. prob from all states in the previous frame to al states
    trans = compute_trans_prob_multichannel(frame, combinations);

    // the index of the final state of the best previous cluster
    int bestGlobalState = R[frame-1]*min_duration + min_duration-1;

    // 
    // Handle state 0 of each model
    //   This state *must* come from the final state of best model
    //   at the previous time.
    // It only makes sense if we have a minimum duration>1
    //
    if(min_duration > 1)
    {
      for(int cl=0; cl<nClass; cl++){
        int curState;
        // we compute the transition prob for that cluster from the previous best
        // to the current computed.
        double clusterTrans = dC[LAST][bestGlobalState] + trans_weight * trans[R[frame-1]][cl];
        //
        // Get the index of the 0th state of current cluster
        //
        curState = cl*min_duration;	
	  
        // add local distance to state's global score
        dC[CURR][curState] = clusterTrans + (double)(lkld[cl]);
        // since we are transitioning now, this frame is the begin time
        tC[CURR][curState] = frame;
      }
    }

    // 
    // Handle the other states of each cluster
    //
    for(int cl=0; cl<nClass;cl++) {
      double forwardTrans;
      int gl_state;		// index into global state array
      int cl_state;		// index of state in a cluster
      double localDist = (double)(lkld[cl]);
      

      //
      // Start with state 1 and end at the second-to-last state
      //   There are no self-loops on the internal states!
      //
      for(cl_state = 1; cl_state < (min_duration-1); cl_state++) 
      {
        gl_state = (cl*min_duration) + cl_state;
        //WARNING: We consider a trans probability so that it counts in each
        // step of the process, regardless of the min. duration.
        //it could be understood as an addendum to the state lkld
        forwardTrans = dC[LAST][gl_state-1] + trans_weight * trans[cl][cl]; 
        dC[CURR][gl_state] = forwardTrans + localDist;
        tC[CURR][gl_state] = tC[LAST][gl_state-1];
      }
      
      //
      // Now do the final state of each cluster
      //  This is handled differently from the internal states because
      //  it has a self-loop prob
      // this changes slightly if min. dur. is 1
      //
      if(min_duration > 1)
      {
        gl_state = (cl*min_duration) + cl_state;
        forwardTrans = dC[LAST][gl_state-1] + trans_weight * trans[cl][cl];
      }
      else
      {
        //I come from the best in the prev run
        gl_state = cl;
        forwardTrans = dC[LAST][bestGlobalState] + trans_weight * trans[R[frame-1]][cl];
      }
      double selfLoopTrans = dC[LAST][gl_state] + trans_weight * trans[cl][cl];
      if (selfLoopTrans >= forwardTrans) {
	dC[CURR][gl_state] = selfLoopTrans + localDist;
	tC[CURR][gl_state] = tC[LAST][gl_state];
      } else {
	dC[CURR][gl_state] = forwardTrans + localDist;
	if(min_duration>1)
        {
          tC[CURR][gl_state] = tC[LAST][gl_state-1];
        }
	else
        {
          //we transition now, we point to here as the start of the cluster
          tC[CURR][gl_state] = frame;
        }
      }
      //test
      if(dC[CURR][gl_state] < -1e30)
      {
        char tmp_string[1024]; sprintf(tmp_string, "WARNING: value dc for state %d is %f\n", gl_state, dC[CURR][gl_state]);
        m_log.print_this(tmp_string, 5);
      }
    }
    
    // 
    // Finished with all states of all models. Now find the cluster
    // that has the highest likelihood in its final state.
    //
    double MaxLkld = -1e30;
    int MaxClass = -1;
    int MaxClassFinSt = -1;
    for(int cl=0; cl<nClass; cl++) {
      int finalClusState = cl*min_duration + min_duration-1;
      double tmplkld = dC[CURR][finalClusState];
      if(tmplkld > MaxLkld) {
	MaxClass = cl;
	MaxLkld = tmplkld;
	MaxClassFinSt = finalClusState;
      }
    }
    
    if (MaxClass == -1) {
      m_log.print_this("ERROR: error finding best class\n", 10);
      exit(1);
    }
    // Record the results
    R[frame] = MaxClass;
    F[frame] = tC[CURR][MaxClassFinSt];

    // copy CURR to LAST
    for(int st=0;st <totStates; st++) {
      dC[LAST][st] = dC[CURR][st];
      tC[LAST][st] = tC[CURR][st];
    }
  } 
  
  //record the best delays on that selected case
  //last one
  int from = F[nSamples-1];
  int to;
  int position;
  for(int i = from; i<nSamples; i++)
  {
    //each of the positions of the winning vector gives a delay position
    position=0;
    for(int channel_count = 0; channel_count<m_numCh; channel_count++)
    {
      //if(channel_count != (*m_config).reference_channel)
      //{
        m_bestClass[channel_count][i] = combinations[R[nSamples-1]][position];
        position++;
      //}
    }
  }
  //others
  while(from>1)
  {
    to = from - 1; //this one included
    from = F[to];
    for(int i=from; i<=to; i++)
    {
      //each of the positions of the winning vector gives a delay position
      position=0;
      for(int channel_count = 0; channel_count<m_numCh; channel_count++)
      {
        //if(channel_count != (*m_config).reference_channel)
        //{
          m_bestClass[channel_count][i] = combinations[R[to]][position];
          position++;
        //}
      }
    }
  }	  
  
  //print out some results
  int bestFinalClus = R[nSamples-1];
  int bestFinalState = bestFinalClus * min_duration + min_duration - 1;
  double viterbi_score = (dC[CURR][bestFinalState]); //I cast back from double to float
  sprintf(tmp_string, "   Segmentation Over with Viterbi score = %f\n", viterbi_score);
  m_log.print_this(tmp_string, 1);
  fflush(stdout);
  
}

/*!
  Computes the transition probability in the mono-channel case

  \param frame segment nunber being computed
  \param channel_count Channel being considered
  \return Transition probabilities matrix between all TDOA N-best pairs in that segment

*/

std::vector<vector<double> > TDOA::compute_trans_prob(int frame, int channel_count)
{
  // the transition prob is defined as abs(curr_delay - prev_delay) in prob. level
  //need to normalize the value, invert it and normalize again, and log it.
  double max_trans;
  vector<vector<double> > trans;
  double LOG_10 = log(10.0);

  for(unsigned int i=0; i<(*m_config).nbest_amount; i++)
  {
    trans.push_back(vector<double>((*m_config).nbest_amount,0));
  }

  //compute the matrix of distances and their total sum
  //total_trans_sum = 0;
  max_trans = -1000;
  for(unsigned int i=0; i<(*m_config).nbest_amount; i++)
  {
    for(unsigned int j=0; j<(*m_config).nbest_amount; j++)
    {
      trans[i][j] = (double)(abs(m_chDelays[channel_count][frame][j] - m_chDelays[channel_count][frame-1][i]));
      if(trans[i][j]>max_trans){max_trans = trans[i][j];}
    }  
  }
  
  //change scale, normalize and log
  /*
  printf("Transition prob:\n");
  for(unsigned int j=0; j<(*m_config).nbest_amount; j++)
  {
    printf(" %d ", m_chDelays[channel_count][frame][j]);
  }  
  printf("\n");
  */
  for(unsigned int i=0; i<(*m_config).nbest_amount; i++)
  {
    //printf("%d: ", m_chDelays[channel_count][frame-1][i]);
    for(unsigned int j=0; j<(*m_config).nbest_amount; j++)
    {
      trans[i][j] = log((max_trans +1 - trans[i][j])/(max_trans+2))/LOG_10;
      //printf(" %f ", trans[i][j]);
    }  
    //printf("\n");
  }
  
  return(trans);
}


/*!
  Computes the transition probability for the multichannel decoding. To do so, it computes the trans. prob. between all allowed combinations
  of channel n-best values. This way we can have better control wether computing all combinations of there are some that are not possible.

  \param frame Number of the segment processed
  \param combinations Vector of the possible combinations
  \return A matrix with all transition probabilities between all possible combinations

*/

std::vector<vector<double> > TDOA::compute_trans_prob_multichannel(int frame, vector<vector<int> > combinations)
{
  // the transition prob is defined as abs(curr_delay - prev_delay) in prob. level
  //need to normalize the value, invert it and normalize again, and log it.
  double max_trans;
  double LOG_10 = log(10.0);

  //where I'll put the output trans
  vector<vector<double> > trans;
  for(unsigned int i=0; i<combinations.size(); i++)
  {
    trans.push_back(vector<double>(combinations.size(),0));
  }

  //where I'll put the local trans
  vector<vector<vector<double> > > tmp_trans;
  vector<vector<double> > tmp_tmp;
  for(unsigned int i=0; i<(*m_config).nbest_amount; i++)
  {
    tmp_tmp.push_back(vector<double>((*m_config).nbest_amount,0));
  }
  for(int i=0; i<(m_numCh-1); i++)
  {
    tmp_trans.push_back(tmp_tmp);
  }  


  //compute the matrix of distances and their total sum
  //total_trans_sum = 0;
  max_trans = -1;
  for(int channel_count = 0; channel_count< m_numCh-1; channel_count++)
  {
    for(unsigned int i=0; i<(*m_config).nbest_amount; i++)
    {
      for(unsigned int j=0; j<(*m_config).nbest_amount; j++)
      {
        tmp_trans[channel_count][i][j] = (double)(abs(m_chDelays[channel_count][frame][j] - m_chDelays[channel_count][frame-1][i]));
        if(tmp_trans[channel_count][i][j]>max_trans){max_trans = tmp_trans[channel_count][i][j];}
      }  
    }
  }

  //change scale, normalize and log
  for(int channel_count = 0; channel_count< m_numCh-1; channel_count++)
  {
    for(unsigned int i=0; i<(*m_config).nbest_amount; i++)
    {
      for(unsigned int j=0; j<(*m_config).nbest_amount; j++)
      {
        tmp_trans[channel_count][i][j] = log((max_trans +1 - tmp_trans[channel_count][i][j])/(max_trans+2))/LOG_10;
      }  
    }
  } 
 
  //now we fill up the output trans prob vector, considering independence between all probs.
  for(unsigned int comb1=0; comb1<combinations.size(); comb1++)
  {
    for(unsigned int comb2=0; comb2<combinations.size(); comb2++)
    {
      //in each position we add up all the trans probabilities
      for(int i=0; i<m_numCh-1; i++)
      {
        trans[comb1][comb2] += tmp_trans[i][combinations[comb1][i]][combinations[comb2][i]];
      }
    }
  }

  return(trans);
}

/*!
  compute the prob. for each element by normalizing to 1 the sum of all, and taking the log. Done for mono-channel case
  
  \param frame Number of th segment computed
  \param channel_count Channel taken into account
  \return probability vector for all n-best values

*/

std::vector<double> TDOA::compute_xcorr_prob(int frame, int channel_count)
{
  double total_lkld_sum=0;
  vector<double> lkld((*m_config).nbest_amount,0);
  double LOG_10 = log(10.0);

  for(unsigned int i=0; i<(*m_config).nbest_amount; i++)
  {
    total_lkld_sum += m_chXcorrValues[channel_count][frame][i];
  }
  for(unsigned int i=0; i<(*m_config).nbest_amount; i++)
  {
    if(m_chXcorrValues[channel_count][frame][i]>0)
    {
      lkld[i] = log((double)(m_chXcorrValues[channel_count][frame][i])/total_lkld_sum)/LOG_10;
    }
    else
    {
      //when the prob is 0
      lkld[i] = -1000;//log(1e-30)/LOG_10;
    }
    //test
    //printf("Total lkld %f, Lkld channel %d: %f\n", total_lkld_sum, i, lkld[i]); 
  }
  
  return(lkld);
}

/*!
  computes the probability of all combinations of multichannel delays given their xcorr
  the possible combinations are passed as a parameter
  
  \param frame Number of th segment computed
  \param combinations list of the n-best combinations that we want to compute the probability for
  \return probability vector for all combinations
 
*/

std::vector<double> TDOA::compute_xcorr_prob_multichannel(int frame, vector<vector<int> > combinations)
{  
 
  vector<double> lkld(combinations.size(),0);
  float LOG_10 = log(10.0);
  int position;

  for(unsigned int i=0; i< combinations.size(); i++)
  {      
    //compute the lkld
    position = 0;
    lkld[i] = 0;
    for(int j=0; j<m_numCh; j++)
    {
      //if(j != (*m_config).reference_channel)
      //{
        if(m_chXcorrValues[j][frame][combinations[i][position]] > 0)
        {
          lkld[i] += log(m_chXcorrValues[j][frame][combinations[i][position]])/LOG_10;
        }
        position++;
      //}	  
    }
  }

  return(lkld);
}

/*!
  creates a vector with all combinations of possible delays each vector has a different 
  combination of the Nbest delay-positions.
  It creates only the possible combinations for the specified block of channels
  for all others it gets the optimum values (1st best). The combinations vector contains the position
  in the TDOA vector for each channel that needs to be used, according to the 1/2-best values computed in
  the mono-channel pass.

  \param Nbest
  \param block_size size of the block of channels for which combinations are created
  \param offset Offset from 0 that the combinations are created
  \param frame segment number to compute the combinations for
  \return the list of combinations computed
*/
std::vector<vector<int> > TDOA::create_delays_positions_block(int Nbest, int block_size, int offset, int frame)
{
  //int possible_combin = (int)(pow((double)(Nbest), (double)(m_numCh-1))); //total amount of combinations
  int possible_combin = (int)(pow((double)(Nbest), (double)(block_size))); //total number of combinations

  //first dimension is the vector number, second dimension are the positions
  //we create a sub-array with the permuted values
  vector<vector<int> > positions;
  for(int i=0; i<possible_combin; i++)
  {
    positions.push_back(vector<int>(block_size,0));
  }

  for(int i=1; i<possible_combin; i++)
  {
    //first we copy the previous value
    for(int j=0; j<block_size-1; j++)
    {
      positions[i][j] = positions[i-1][j];
    }
    positions[i][block_size-1] = positions[i-1][block_size-1] + 1; //previous one +1 in the last value

    //we work from backwards to forward rearranging all values
    for(int j=block_size-1; j>0; j--)
    {
      if(positions[i][j] >= Nbest)
      {
        positions[i][j] = 0;
        positions[i][j-1]++;
      }
    }
  }

  //now we fill-out the real vector with extended values according to the best value
  //vector<int> best_positions_vector((m_numCh-1),0);
  vector<int> best_positions_vector(m_numCh, 0);
  int position=0;
  for(int channel_count = 0; channel_count<m_numCh; channel_count++)
  {
    //if(channel_count != (*m_config).reference_channel)
    //{
      best_positions_vector[position] = m_bestClass[channel_count][frame];
      position++;
    //}
  }
  
  vector<vector<int> > positions_final;
  for(int i=0; i<possible_combin; i++)
  {
    positions_final.push_back(vector<int>(best_positions_vector));
  }

  for(int i=0; i<possible_combin; i++)
  {
    //permuted values
    for(int j=0; j<block_size; j++)
    {
      positions_final[i][offset+j] = positions[i][j];
    }
  }

  return(positions_final);
}

/*******************************************************************************
* Continuity filtering
********************************************************************************/

/*!
  Applies a simple acoustic continuity filter to the signal, based on the n-best TDOA
  values computed and the previous value selected.
  
  \return Nothing
*/

void TDOA::continuity_filter()
{

  //tells us wether the next delay is found to be continuous with the previous
  int CONTINUOUS_DELAY ;
  
  for(int channel_count=0; channel_count < m_numCh; channel_count++)
  {
    //the current algorithms with varying reference channel need this different
    //if(channel_count != (*m_config).reference_channel)
    //{
      //go through all computed delays
    for(int count=0; count<(int)(m_chDelays[channel_count].size()); count++)
    {
      //check whether that channel delay is the reference
      if(m_chXcorrValues[channel_count][count][0] != 1)
      {

        CONTINUOUS_DELAY = 1;
	      
        if(count>0)
        {
          //we only do anything when we are not at the first frame
		  
          //delay continuity
          if((m_chDelays[channel_count][count-1][0] + (*m_config).delay_variance) < m_chDelays[channel_count][count][0] || 
             (m_chDelays[channel_count][count-1][0] - (*m_config).delay_variance) > m_chDelays[channel_count][count][0])
          {
            CONTINUOUS_DELAY = 0;
            //the first one doesn't follow continuity, we see if there is any that does
            for(unsigned int i=1; i < (*m_config).nbest_amount; i++)
            {
              if((m_chDelays[channel_count][count-1][0] + (*m_config).delay_variance) > m_chDelays[channel_count][count][i] && 
                 (m_chDelays[channel_count][count-1][0] - (*m_config).delay_variance) < m_chDelays[channel_count][count][i])
              {
                //enable the flag
                m_delayFilters[channel_count][count] += F_CONTIN;

                //we found it, we rearrange the order of the delays and move on
                //we move all down from the position "i" til the first
                int tmp_delay = m_chDelays[channel_count][count][i];
                float tmp_xcorr_value = m_chXcorrValues[channel_count][count][i];
                for(int j=i; j>0; j--)
                {
                  m_chDelays[channel_count][count][j] = m_chDelays[channel_count][count][j-1];
                  m_chXcorrValues[channel_count][count][j] = m_chXcorrValues[channel_count][count][j-1];				  
                }
                m_chDelays[channel_count][count][0] = tmp_delay;
                m_chXcorrValues[channel_count][count][0] = tmp_xcorr_value;
                CONTINUOUS_DELAY = 1;
                break;
              }
            }
          }

        }	      
      }
    }//end of the delay frames loop
  }//channels loop
}

/*!
  compute the n-best GCC-PHAT values. It is an obsolete function now.
  \param counter Determines the starting point of the analysis window
  \param chanData channel data to compute the channel's xcorr from
  \return Nothing
*/

void TDOA::compute_channels_xcorr(long counter, float** chanData)
{
  //holders for the delays info for 1 channel
  vector<int> delays((*m_config).nbest_amount,0);
  vector<float> xcorr_values((*m_config).nbest_amount,0);

  for(int channel_count=0; channel_count< m_numCh; channel_count++)
  {
    if(channel_count != (*m_config).reference_channel)
    {
      xcorrelation_FFTReal(&delays[0], &xcorr_values[0], (*m_config).nbest_amount, chanData[channel_count], chanData[(*m_config).reference_channel], (*m_config).marginFrames, (*m_config).windowFrames, m_fftWindowFrames);
    }
      
    else
    {
      //for the reference channel, the delay is set to 0 and we calculate the power (xcorr)
      delays.assign((*m_config).nbest_amount, 0);
      xcorr_values.assign((*m_config).nbest_amount, 0);
      xcorr_values[0] = 1;
    }

    //we insert the values into the general structure
    m_chDelays[channel_count][counter] = delays;
    m_chXcorrValues[channel_count][counter] = xcorr_values;

    //we print out the values
    char tmp_string[1024]; sprintf(tmp_string," %d (%.3f)\t", m_chDelays[channel_count][counter][0], m_chXcorrValues[channel_count][counter][0]);
    m_log.print_this(tmp_string, 0);
  }	 
}


/*!
  Prints out both the delays and xcorr values finally used on the channels
  
  \\param SampleRateInMs Sample Rate in Milliseconds
  \return Nothing
*/

void TDOA::print_delays(int sampleRateInMs)
{
  //prints the delays and xcorr values
  //looks at the m_delayFilters to print information about the kind of filterings done
  long frame;

  for(int count=0; count<(int)(m_chDelays[0].size()); count++)
  {
    //print the line header for output delays
    frame = long((*m_config).rate * sampleRateInMs * count);
    fprintf((*m_fileInOut).delfd,"%ld ",frame);
    fprintf((*m_fileInOut).delfd2,"%ld ",frame);
      
    for(int channel_count=0; channel_count < m_numCh; channel_count++)
    {
      //go through all computed delays
      fprintf((*m_fileInOut).delfd2,"%d %f ", m_chDelays[channel_count][count][0], m_chXcorrValues[channel_count][count][0]);
      fprintf((*m_fileInOut).delfd,"%d %f ", m_chDelays[channel_count][count][1], m_chXcorrValues[channel_count][count][1]);	 
    }
    fprintf((*m_fileInOut).delfd, "\n");
    fprintf((*m_fileInOut).delfd2, "\n");  

  }  
}

/*!
  does the GCC-PHAT optimizing for processing speed. If computes all channels keeping the acoustic data and some FFT's from one to the other.

  \param counter Segment number being processed
  \param chanData Channels data
  \return the 1best delays obtained
*/

vector<int> TDOA::xcorrelation_FFTReal_full(long counter, vector<vector<float> > chanData)
{

    int complexSize = m_fftWindowFrames; //we use an actual FFT of 2*m_fftWindowFrames
    //we first compute the FFT for the input channels, all at once so that the FFT of the reference
    //channel is computed just once

    //instantiate the FFT
    ffft::FFTReal <float> fft_object (2*m_fftWindowFrames);

    //1) create the containers and assign the input data
    float** channels_fft = new float*[m_numCh];
    float** channels_in = new float*[m_numCh];

    vector<float> hamm_val((*m_config).windowFrames);
    for(int i=0; i<(*m_config).windowFrames; i++)
        hamm_val[i] = 0.54 - 0.46*cos(6.283185307*i/((*m_config).windowFrames-1));

    for(int chan_count=0; chan_count<m_numCh; chan_count++)
    {
        channels_in[chan_count] = new float[2 * m_fftWindowFrames];

        //we need to do a padding of 0'os at the end
        //we multiply by a hamming window
        for(int i=0; i < (*m_config).windowFrames; i++)
            channels_in[chan_count][i] = (float)(chanData[chan_count][i]) * hamm_val[i];
        for(int i = (*m_config).windowFrames; i < (2*m_fftWindowFrames); i++)
            channels_in[chan_count][i] = 0.00;

        //perform the FFT
        channels_fft[chan_count] = new float[2 * m_fftWindowFrames];
        fft_object.do_fft (channels_fft[chan_count], channels_in[chan_count]);
    }

  //We now have the FFT values for each channel
  // we compute the GCC-PHAT and search the maxima, for the desired combination of channel pairs
  float *result;
  result = new float[2*m_fftWindowFrames];
  float *xcorr_value;
  xcorr_value = new float[2*(*m_config).marginFrames+1];

  float* tmpData = new float[2*m_fftWindowFrames];
  float abs_value;
  int refChannel = (*m_config).reference_channel;
  for(int channel_count=0; channel_count< m_numCh; channel_count++)
  {
    if(channel_count != (*m_config).reference_channel)
    {
        //here we do the GCC-Phat computation (only half needed)
        //we do it directly over the output of the FFTReal, taking into account the way it structures the data
        for(int i=0; i <= complexSize; i++)
        {
            //do the multiplication: (a+jb)(c+jd)* = (ac+bd) - j(ad+bc)
            if(i==0 || i==complexSize)
            {
                tmpData[i] = channels_fft[channel_count][i] * channels_fft[refChannel][i];
                abs_value = sqrt(tmpData[i]*tmpData[i]);
            }
            else
            {
                tmpData[i]             =  channels_fft[channel_count][i] * channels_fft[refChannel][i] + channels_fft[channel_count][complexSize+i] * channels_fft[refChannel][complexSize+i];
                tmpData[complexSize+i] =  channels_fft[channel_count][complexSize+i] * channels_fft[refChannel][i] - channels_fft[channel_count][i] * channels_fft[refChannel][complexSize+i];
                abs_value = sqrt(tmpData[i]*tmpData[i] + tmpData[complexSize+i]*tmpData[complexSize+i]);
            }

            //divide by the abs value
            if(abs_value == 0)
            {
                //avoid division by 0;
                abs_value = 1;
            }

            //we divide and normalize the extra scaling
            tmpData[i] /= (abs_value * 2 * m_fftWindowFrames);
            if(i != complexSize && i != 0)
                tmpData[complexSize+i] /= (abs_value * 2 * m_fftWindowFrames);
        }
        //do the inverse fft of the result
        fft_object.do_ifft(tmpData, result);

        //we assign the result to the xcorr_value variable. Result is a "real" vector
        //We need to do a shift of the data, as we want positive and negative parts
        //positive part:
        for(int i=0; i<(*m_config).marginFrames; i++)
        {
            //result[i] /= 2*m_fftWindowFrames;
            xcorr_value[i+((*m_config).marginFrames+1)] = (float)(result[i]);
        }
        //negative part
        for(int i=2*m_fftWindowFrames-((*m_config).marginFrames+1); i<2*m_fftWindowFrames; i++)
        {
            //result[i] /= 2*m_fftWindowFrames;
            xcorr_value[i-(2*m_fftWindowFrames-((*m_config).marginFrames+1))] = (float)(result[i]);
        }

        //in here the normal stuff continues: find the N-best values and set then into the delay vectors
        find_nbest_maximums(&m_chDelays[channel_count][counter][0], &m_chXcorrValues[channel_count][counter][0], (*m_config).nbest_amount, xcorr_value, (*m_config).marginFrames, m_fftWindowFrames);
    }
    else //when we have the reference channel
    {
        //for the reference channel, the delay is set to 0 and we calculate the power (xcorr)
        m_chDelays[channel_count][counter].assign((*m_config).nbest_amount, 0);
        m_chXcorrValues[channel_count][counter].assign((*m_config).nbest_amount, 0);
        m_chXcorrValues[channel_count][counter][0] = 1;
    }
  }
  //clean up
  delete [] tmpData;


  //Independently of the reference channel computation, always ensure continuity in the reference channel versus its previous value
  //Done this way so that reference channel computation can be flexibilized
  vector<int> delays_out;
  //int tmpDelays = m_chDelays[0][counter][0]; //save the delays for channel 0
  //m_chDelays[0][counter].assign((*m_config).nbest_amount, 0);

  if((*m_config).COMPUTE_REFERENCE == 2)
  {
        //independently of the channel chosen as reference, we move the delays so that channel 0 has 0 delay
        int delay = m_chDelays[0][counter][0]; //the first delay of channel 0
        printf("reference %d Delay %d\n", (*m_config).reference_channel, delay);
        for(int i=0; i<m_numCh; i++)
            for(unsigned int nbest=0; nbest<(*m_config).nbest_amount; nbest++)
        m_chDelays[i][counter][nbest] -= delay;

        /*
    //We compute all possible alignments between channels and compute the sum of resulting differences. We then select that which minimizes the sum
    int best_channel = (*m_config).reference_channel; //by default the reference channel
    int distances_sum = 0;
    int best_distances = 1000000;
    if(counter > 0)
    {
      //find which gap reduces it the most
      int chan_gap;
      for(int i=0; i<m_numCh; i++)
      {
        chan_gap = m_chDelays[i][counter][0] - m_chDelays[i][counter-1][0];
        //compute current distances
        distances_sum = 0;
        for(int channel_count=0; channel_count< m_numCh; channel_count++)
        {
          distances_sum += abs(m_chDelays[channel_count][counter][0] - m_chDelays[channel_count][counter-1][0] - chan_gap);
        }
        //printf("Distances sum channel %d: %d\n", i, distances_sum);

        if(distances_sum < best_distances)
        {
          best_channel = i;
          best_distances = distances_sum;
        }
      }
      //printf("Synching delays with channel %d\n", best_channel);
    }

    for(int channel_count=0; channel_count< m_numCh; channel_count++)
    {
      if(counter > 0)
      {
        for(unsigned int numBest = 0; numBest < (*m_config).nbest_amount; numBest++)
        {
          m_chDelays[channel_count][counter][numBest] -= m_chDelays[best_channel][counter][0] - m_chDelays[best_channel][counter-1][0];
        }
      }
    }
        */
  }
/*
  for(int channel_count=0; channel_count< m_numCh; channel_count++)
  {
    if(counter > 0)
    {
      for(int numBest = 0; numBest < (*m_config).nbest_amount; numBest++)
      {
        //add the previous reference delay
        m_chDelays[channel_count][counter][numBest] += m_chDelays[(*m_config).reference_channel][counter-1][0];
      }
    }
*/

    //print out the values computed
  for(int channel_count=0; channel_count< m_numCh; channel_count++)
  {
    //we print out the values
    char tmp_string[1024]; sprintf(tmp_string, " %d (%.3f)\t", m_chDelays[channel_count][counter][0], m_chXcorrValues[channel_count][counter][0]);
    m_log.print_this(tmp_string, 0);
    delays_out.push_back(m_chDelays[channel_count][counter][0]);
  }

  //now we can free the rest of variables
  //we free and destroy the used variables
  delete [] result;
  delete [] xcorr_value;
  for(int i=0; i<m_numCh; i++)
  {
      delete [] channels_fft[i];
      delete [] channels_in[i];
  }


  return(delays_out);
}



/*!
  Given a xcorrelation function, finds the n-best maximum values and their locations
  \param delays
  \param delays TDOA values vector
  \param values TDOA values xcorr vector
  \param amount_max Number of N-best TDOA values returned
  \param xcorr_value Vector containing the xcorr_values to find the maximum points of
  \param margin Number of frames around frame 0 that delays are looked for (sets the maximum possible delay to be found)
  \param fftwindow_size Number of frames of the FFT
  \return Sum of the xcorr_values across all data, indicative of the quality of the comparison.
*/

float TDOA::find_nbest_maximums(int *delays, float *values, int amount_max, float *xcorr_value, int margin, int fftwindow_size)
{
    (void)fftwindow_size; //this is unused

  //maximums masking value
  int masking = 5;
  
  //we find the maximums of the function and return amount_maximums
  vector<float> max_values;
  vector<int> max_pos;
  
  //we add all the xcorr values to return the overall area
  float xcorr_sum;
  xcorr_sum = xcorr_value[0];
  //finding the maximums
  //I don't consider the delay=0 a maximum, as it has a discontinuity always.
  for(int i=1; i<(2*margin)-1; i++)
  {
    xcorr_sum += xcorr_value[i];
    if(xcorr_value[i]>xcorr_value[i-1] && xcorr_value[i]>xcorr_value[i+1])
    {
      //we have a maximum
      max_values.push_back(xcorr_value[i]);
      max_pos.push_back(i);
    }
  }
  
  //check that there is at least one maximum. If not, we set the delay to 0 with very low value
  if(max_values.size() == 0)
  {
    max_values.push_back(0.0001);
    max_pos.push_back(margin);
    m_log.print_this("WARNING!!! no maximum points found in crosscorrelation!!! does your signal contain digital zeros?\n", 5);
  }
  
  //we find the amount_max maximums from the maximums list
  float max_value=-1;
  int max_idx=-1;
  int max_idx_here=-1;
  int good_max;
  for(int i=0; i<amount_max; i++)
  {
    max_value = -1;
    for(int count=0; count<(int)(max_values.size()); count++)
    {
      if(max_values[count] > max_value)
      {
        max_value = max_values[count];//the value of the maximum
        max_idx = max_pos[count];//the delay of the maximum
        max_idx_here = count;//the position in the vector of maximums
      }
    }
      
    //we check that no other preselected maximum is too close to this one
    //the masking margin around a maximum is set empirically to 5
    if(max_value != -1)
    {
      good_max=1;
      for(int j=0; j<i; j++)
      {
        if(max_idx-masking<(delays[j]+margin+1) && max_idx+masking>(delays[j]+margin+1))
        {
          //we already have a maximum around there
          good_max = 0;
          max_values[max_idx_here] = -1;//we don't want this one to be looked again
        }
      }
      if(good_max == 0)
      {
        i--; //we need to find another one	  
      }
      else
      {
	      
        delays[i] = max_idx - margin - 1; //delay relative to +-margin
        values[i] = max_value;
        max_values[max_idx_here] = -1;
      }
    }
    else
    {
      //we didn't find any more maximums, we repeat the first delay and value (it's a way to put null)
      delays[i]=delays[0];
      values[i] = values[0];
    }
  }
  if(max_idx == -1 || max_idx_here == -1)
  {
    m_log.print_this("ERROR: Maximum index not found\n", 10);
    exit(1);
  }

  return(xcorr_sum);
}


/*!
  we do the xcorrelation using the FFTReal library functions. Only the GCC-Phat is implemented at this stage

  \param delays TDOA values vector
  \param values TDOA values xcorr vector
  \param amount_max Number of N-best TDOA values returned
  \param chan_data Channel vector data
  \param ref_data Reference channel vector data
  \param margin Number of frames around frame 0 that delays are looked for (sets the maximum possible delay to be found)
  \param window Size (in frames) of the analysis window
  \param fftwindow_size Number of frames of the FFT
 */

float TDOA::xcorrelation_FFTReal(int *delays, float *values, int amount_max, float *chan_data, float *ref_data, int margin, int window, int fftwindow_size)
{

  int i;
  //int freq_filter;
  float *xcorr_value;
  float hamm_val;
  xcorr_value = new float[2*margin+1];
  int complexSize = fftwindow_size; //we use an actual FFT of 2*fftwindow_size

  //create the fft object (it would be faster if created only once)
  //TODO: create it whenever I have the window size known. Consider using the fixed length implementation
  ffft::FFTReal <float> fft_object (2*fftwindow_size);

  //reserve the data
  //TODO: check whether I can do this outside, as it takes time and I always reserve the same space
  float *channel_in, *refer_in, *result, *channelFFT, *refFFT;
  result = new float[2*fftwindow_size];
  channel_in = new float[2*fftwindow_size];
  refer_in = new float[2*fftwindow_size];
  channelFFT = new float[2*fftwindow_size];
  refFFT = new float[2*fftwindow_size];

  //we assign the input data:
  //we assign the data to channel and refer
  //we need to do a padding of 0'os at the end
  //we multiply by a hamming window
  for(i=0; i<(window); i++)
  {
    hamm_val = 0.54 - 0.46*cos(6.283185307*i/(window-1));
    channel_in[i] = (float)(chan_data[i]) * hamm_val;
    refer_in[i] = (float)(ref_data[i]) * hamm_val;
  }
  for(i=window; i<(2*fftwindow_size); i++)
  {
    channel_in[i] = 0;
    refer_in[i] = 0;
  }

  //execute the FFT

  //TODO: put these outside of the function
  fft_object.do_fft (channelFFT, channel_in);
  fft_object.do_fft (refFFT, refer_in);

  //here we do the GCC-Phat computation (only half needed)
  //we do it directly over the output of the FFTReal, taking into account the way it structures the data
  float* tmpData = new float[2*fftwindow_size];
  float abs_value;
  for(i=0; i <= complexSize; i++)
  {
      //do the multiplication: (a+jb)(c+jd)* = (ac+bd) - j(ad+bc)
      if(i==0 || i==complexSize)
      {
          tmpData[i] = channelFFT[i] * refFFT[i];
          abs_value = sqrt(tmpData[i]*tmpData[i]);
      }
      else
      {
          tmpData[i] = channelFFT[i] * refFFT[i] + channelFFT[complexSize+i] * refFFT[complexSize+i];
          tmpData[complexSize+i] =  channelFFT[complexSize+i] * refFFT[i] - channelFFT[i] * refFFT[complexSize+i];
          abs_value = sqrt(tmpData[i]*tmpData[i] + tmpData[complexSize+i]*tmpData[complexSize+i]);
      }

      //divide by the abs value
      if(abs_value == 0)
      {
          //avoid division by 0;
          abs_value = 1;
      }

      //we divide and normalize the extra scaling
      tmpData[i] /= (abs_value * 2 * fftwindow_size);
      if(i != complexSize && i != 0)
          tmpData[complexSize+i] /= (abs_value * 2 * fftwindow_size);
  }

  //do the inverse fft of the result
  fft_object.do_ifft(tmpData, result);
  delete [] tmpData;

  //we assign the result to the xcorr_value variable. Result is a "real" vector
  //We need to do a shift of the data, as we want positive and negative parts
  //positive part:
  for(i=0; i<margin; i++)
  {
    //result[i] /= 2*fftwindow_size;
    xcorr_value[i+(margin+1)] = (float)(result[i]);
  }
  //negative part
  for(i=2*fftwindow_size-(margin+1); i<2*fftwindow_size; i++)
  {
    //result[i] /= 2*fftwindow_size;
    xcorr_value[i-(2*fftwindow_size-(margin+1))] = (float)(result[i]);
  }

  //in here the normal stuff continues
  float xcorr_sum;
  xcorr_sum = find_nbest_maximums(&delays[0], &values[0], amount_max, xcorr_value, margin, fftwindow_size);

  //we free and destroy the used variables
  delete [] channel_in;
  delete [] refer_in;
  delete [] result;
  delete [] xcorr_value;
  delete [] channelFFT;
  delete [] refFFT;

  //return
  return(xcorr_sum);

}


/*!
  writes out the finalized optimum delays.
  \return Nothing
*/

void TDOA::Optimum_delays_write_out()
{
  //once filtered I copy all of them into the final structures
  //obtaining the best delay via the m_bestClass vector
  for(int channel_count=0; channel_count < m_numCh; channel_count++)
  {
    m_finalDelays[channel_count].clear();
    m_finalXcorrValues[channel_count].clear();
    for(int count=0; count < (int)(m_chDelays[channel_count].size()); count++)
    {
      //we write down the first and second best delays
      m_finalDelays[channel_count].push_back(vector<int>(2));
      m_finalXcorrValues[channel_count].push_back(vector<float>(2));
	 	 
      if(m_bestClass[channel_count][count] == 0)
      {
        m_finalDelays[channel_count][count][0] = m_chDelays[channel_count][count][0];
        m_finalXcorrValues[channel_count][count][0] = m_chXcorrValues[channel_count][count][0];

        m_finalDelays[channel_count][count][1] = m_chDelays[channel_count][count][1];
        m_finalXcorrValues[channel_count][count][1] = m_chXcorrValues[channel_count][count][1];	      
      }
      else
      {
        m_finalDelays[channel_count][count][0] = m_chDelays[channel_count][count][1];
        m_finalXcorrValues[channel_count][count][0] = m_chXcorrValues[channel_count][count][1];

        m_finalDelays[channel_count][count][1] = m_chDelays[channel_count][count][0];
        m_finalXcorrValues[channel_count][count][1] = m_chXcorrValues[channel_count][count][0];	      	      
      }
    }
  }  

  //write to a file the 1st and 2nd best delays found
  m_log.print_this("Printing delays to file\n", 1);
  (*m_fileInOut).delays_to_file(m_finalDelays, m_finalXcorrValues, m_UEMGap);
}

/*!
  returns the 1st best channel delays for a given delays position. It does not check whether these delays have been computed or nor
  \param Nbest number of delays desired
  \return found best channel delays
*/
vector<int> TDOA::get_bestChDelays(int numDelay)
{
  vector<int> tmp_delays;
  for(int i=0; i<m_numCh; i++)
  {
    tmp_delays.push_back(m_chDelays[i][numDelay][0]); 
  }

  return(tmp_delays);
}

