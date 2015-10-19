// Copyright 2006-2014 Xavier Anguera and ICSI.
// Distributed under the ICSI Software License, 
// (See accompanying file license.txt)
// For more information and documentation visit
// www.icsi.berkeley.edu/~xanguera/beamformit


#include <string>
#include <cstring>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <vector>
#include <set>
#include <algorithm>
#include <fstream>

#include "sndfile.h"
//#include "boost/program_options.hpp"
#include "global.h"
#include "delaysum.h"
#include "support.h"
#include "tdoa.h"
#include "fileinout.h"
#include "parse_options.h"

/*! \mainpage BeamformIt documentation page
  \include documentation.txt
*/

//namespace po = boost::program_options;
using namespace std;


int main(int argc, char *argv[]) {

    /*! BeamformIt is composed of this main() routine, 1 class and 2 structures of data
    - DelaySum class: Contains all internal variables and routines necessary to perform the
      filter&sum processing. An instance called "ds" is declared and used.
    - FileInOut structure: Contains the pointers to all file descriptors needed for data
      input and output. It is passed to functions that need to deal with such files.
      The instance used is "fileInOut".
      - Configuration structure: Contains all configuration parameters determined by the user.
        It is instantiated as "config"
    */

	//define the config structure
    Configuration config;
    //define the work_files structure
    FileInOut fileInOut(&config);
    //TDOA computation and processing
    TDOA tdoa(&fileInOut, &config);
    //define the delay_sum instance
    DelaySum ds(&config, &tdoa, &fileInOut);
    //log and Info class
    LogInfo log;

    //////////////////////////////////////////////////////////////////////////////////////////
    //////           SYSTEM INITIALIZATION
    //////////////////////////////////////////////////////////////////////////////////////////

    ///////////deal with the input parameters////////////////
    parse_options poptions(&config);
    poptions.parse_params(argc, argv);

    //po::variables_map vm;
    //ds.input_parameters(argc, argv, vm);

    if(config.print_features)
        log.print_features(config);

    /////////create the necessary directories////////////////////
    log.print_this("Creating necessary directories\n",2);
    //sprintf(config.RESULTPATH, "%s/%s/",config.RESULTPATH.c_str(), config.SHOWNAME.c_str());

    //make sure the output directory string finishes with "/"
    if(config.RESULTPATH.at( strlen(config.RESULTPATH.c_str()) - 1 ) != '/')
      config.RESULTPATH.insert(strlen(config.RESULTPATH.c_str()), "/");
      //config.RESULTPATH.push_back('/');
      //config.RESULTPATH.append("/");
    mkdir(config.RESULTPATH.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);

    /////////////open input/output files//////////////////////
    log.print_this("Opening the input channels\n",2);
    int numChannels = fileInOut.Open_Input_Channels();
    printf("NumChannels is: %d, sampling rate is %ld\n", numChannels, fileInOut.get_sampleRate());fflush(stdout);
    log.print_this("Opening output channels\n",2);
    fileInOut.Open_Output_Channels();

    ///////////set show-specific parameters///////////////////
    ///I can initialize everyting in here, after retrieving all the information from
    /// the config file and the input files
    ds.Set_Show_Parameters();

    //////////if we are using a UEM file, we retrieve the information and set the appropriate values/////////////
    if(config.DO_USE_UEM_FILE)
    {
        char tmp_data[512]; sprintf(tmp_data,"UEM file: %s\n",config.UEMFILE.c_str());
        log.print_this(tmp_data, 2);

        ds.SetUEMlimits();
    }
    else
    {
        ds.set_UEMGap(0);
    }
    //tdoa.set_UEMGap(ds.get_UEMGap());

    //////////////Special test function: write out the residual between two files/////////////
    if(config.DO_OUTPUT_RESIDUAL)
    {
        log.print_this("Computing the residual between two signals\n",2);
        ds.computeChannelsDifference();
        fileInOut.Close_Channels();
        exit(0);
    }

    //////////we compute the skew between signals ////////////
    if(config.DO_COMPUTE_SKEW)
    {
        log.print_this("Computing the skew between input signals",2);
        ds.compute_skew();
    }

    ////////// Initialize the structures to keep all the info, needs to be in the order delaysum + tdoa////////
    log.print_this("Initializing...",2);
    ds.init();
    tdoa.init(numChannels, ds.get_UEMGap());
    log.print_this("finished\n",2);

    ///////////select the reference channel////////////
    ///we select it unless it has been set by configuration
    log.print_this("Setting the reference channel\n",2);
    //compute the Xcorr's and select the best one for the reference
    char tmp_data[512];
    switch(config.COMPUTE_REFERENCE)
    {
    case 0:
        //reference channel selected by default to 0 or to the value determined by the config file
        //check first that the selected reference channel is possible. Exit if not
        if(numChannels <= config.reference_channel)
        {
            sprintf(tmp_data, "ERROR: Number of channels (%d) is smaller than the requested reference (%d), reference channel starts at value 0\n", numChannels, config.reference_channel);
            log.print_this(tmp_data, 10);
            exit(-1);
        }
        sprintf(tmp_data, "Selected channel %d as the reference channel\n", config.reference_channel);
        log.print_this(tmp_data, 2);
        break;
    case 1:
        //reference channel selected at start by xcorr comparison and set for the rest of the recording
        config.reference_channel = ds.find_best_channel_xcorr();
        sprintf(tmp_data, "Selected channel %d as the reference channel\n", config.reference_channel);
        log.print_this(tmp_data, 2);
        break;
    case 2:
        //reference cannel set dynamically thoughout the meting
        config.reference_channel = ds.find_best_channel_xcorr_adapt();
        sprintf(tmp_data, "Dynamic channel selection using xcorr\n");
        log.print_this(tmp_data, 2);
        break;
    default:
        log.print_this("WARNING: Reference channel selection method not selected properly, continuing with default value\n",0);
        break;
    }


    /////////We calculate the scaling factor over the uem chunk (if it exists)////////////
    /// The scaling factor ensures that the dynamic range of the signal fits
    /// optimally into the allocated bytes reserved for it at output
    log.print_this("Computing the scalling factor for the output signal\n",2);
    ds.compute_scalling_factor();
    log.print_this("Finished computation of the scalling factor\n",2);

    ////////////////select the initial channel weights//////////////
    log.print_this("Setting initial Channel Weights...",2);
    ds.Set_Channels_Weights();
    log.print_this("finished",2);


    //////////////////////////////////////////////////////////////////////////////////////////
    //////           TDOA PROCESSING
    //////////////////////////////////////////////////////////////////////////////////////////

    long counter=0;
    long start_frame;

    //char tmp_data[512];
    sprintf(tmp_data, "Frames: %ld\n",ds.get_frames());
    log.print_this(tmp_data, 2);
    sprintf(tmp_data, "Rate %d, window %d frames %ld samplerate %ld\n", config.rate, config.window, ds.get_frames(), ds.get_sampleRate());
    log.print_this(tmp_data, 2);

    /// We compute the N-best TDOA values for each analysis window across the whole recording
    float percent_printed = 0; //percentage printed last
    while(counter < tdoa.get_totalNumDelays())
    {
        //we set the start point for the data
        start_frame = ds.get_start_frame(counter);
        if(start_frame == -1)
        {
            log.print_this("Error: Out of the file\n", 10);
        }

        //print the start time
        float start_time = (float)(start_frame)/ds.get_sampleRate();
        float percent_done = (float)(counter*100)/tdoa.get_totalNumDelays();

        //in order not to clutter the output, we print roughly 100 percentage lines
        if((percent_done - percent_printed) >= 1.0)
        {
            printf("Processing %.3f s. (%.3f percent)\r", start_time, percent_done);fflush(stdout);            
            percent_printed += 1.0;
        }

        //compute the TDOA values
        ds.computeTDOAValues(start_frame, counter);
        counter++;
    }
    log.print_this("\nFinished computing delays\n",2);




    //////////////////////////////////////////////////////////////////////////////////////////
    //////           TDOA VALUES SELECTION
    //////////////////////////////////////////////////////////////////////////////////////////

    //for now we don't reuse the weights computed earlier, we recompute them for each frame
    ds.Set_Channels_Weights();

    /////////// compute the optimum delays from the N-best computed values////////////////
    if(config.DO_OPTIMUM_DELAYS)
    {
        log.print_this("Computing the optimum delays\n",2);

        ///applies filters to the delays
        tdoa.Optimum_delays_filtering();

        ///Selects the best TDOA values
        switch(config.DO_ACOUSTIC_MODELLING)
        {
        case 1:
            //double-pass Viterbi
            tdoa.Optimum_delays_acoustic_modelling();
            break;
        case 2:
            //simple continuity function (much faster)
            tdoa.continuity_filter();
            break;
        default:
            log.print_this("Not applying any continuity filter\n",2);
            break;
        }

        log.print_this("Finished computing the optimum delays\n",2);
    }


    //////////////////////////////////////////////////////////////////////////////////////////
    //////           WRITE OUT RESULTS
    //////////////////////////////////////////////////////////////////////////////////////////


    ///writes out the determined optimum delays
    tdoa.Optimum_delays_write_out();

    ////////// do the sum using the optimum delays found//////////////////
    if(config.DO_WRITE_SPH_FILES)
    {
        log.print_this("Summing up all the channels\n",2);
        ds.Channels_Sum(0); //channel 1
        if(config.output_second_best_sph)
        {
            ds.Channels_Sum(1); //channel 2
        }
        log.print_this("Finished summing up all the channels\n",2);
    }

    ///////////// Optionally, compute a quality measure on the beamformed signal //////////////
    if(config.DO_SIGNAL_QUALITY)
        ds.signalQuality();


    //////////////////////////////////////////////////////////////////////////////////////////
    //////           WRAP UP AND GO HOME
    //////////////////////////////////////////////////////////////////////////////////////////


    /////////////closing the file descriptors/////////////////
    log.print_this("Closing all used channels\n",2);
    fileInOut.Close_Channels();
    log.print_this("Finished closing all the channels\n",2);

}
