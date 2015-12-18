# BeamformIt
#### acoustic beamforming tool


**BeamformIt** is an acoustic beamforming tool that accepts a variable amount of input channels and computes an output via a filter&sum beamforming technique. It makes almost no assumptions on the input data (e.g. number of channels, topology, locations, individual channel audio quality, ...).

BeamformIt was originally implemented by [Xavier Anguera](http://www.xavieranguera.com) at ICSI for participation to the NIST RT05s Meetings evaluation to deal with the different number of microphone channels available in a meeting room. BeamformIt was then rewritten and improved for the RT06s evaluation and finally readjusted and documented for public release.

BeamformIt was initially focused towards processing the data used in the RT evaluations but it can now process all sorts of data. As of version 3.5 an effort has been made to eliminate almost all external library dependencies and a script has been created for "casual" beamforming users to have an easy-to-use tool for their acoustic beamforming needs.

Prior to its release into Github, BeamformIt was released via [this website](http://www.xavieranguera.com/beamformit/) with a versioning system. As of version 3.51 all releases will be done through Github. Following [Kaldi](https://github.com/kaldi-asr/kaldi)'s philosofy, we will not set new versions (but will try to keep a log of key changes [below](#changelog)).

If you use the software for research I would very much appreciate if you could cite my work, you can use any of the following citations:
"Acoustic beamforming for speaker diarization of meetings", Xavier Anguera, Chuck Wooters and Javier Hernando, IEEE Transactions on Audio, Speech and Language Processing, September 2007, volume 15, number 7, pp.2011-2023.
"Robust Speaker Diarization for Meetings", Xavier Anguera, PhD Thesis, UPC Barcelona, 2006.

Index
------

1. [Compiling the code](#compilation)
2. [Running the tool](#execution)
  - [Simple execution](#execute_simple)
  - [(More) complex execution](#execute_complex)
3. [Output files](#output)
4. [How to cite](#cite)
5. [Change Log](#changelog)

<a name="compilation"></a>
Compiling the code
------

As of version 3.5 there is only ONE external library required by BeamformIt:
  - [Libsndfile](http://www.mega-nerd.com/libsndfile/): used to input-output data from the audio files, letting the software deal with the sound files as standard files

In Ubuntu, you can install libsndfile `sudo apt-get install libsndfile1-dev`

Additionally, doxygen can be used to compile the documentation of the source-code

To compile the code, first of all, you need to make sure that the Makefile is pointing to the right directories for the sndfile library (or that it has been installed in the system).
The program uses cmake to compile the code. This means we can usually compile with:
`cmake .`
`make`

Note: if sndfile is not installed in the system (and you do not want to do it) and the compilation complains you can repeate the cmake command as follows:
`cmake -DLIBSND_INSTALL_DIR=/libsnd/install/dir .`

You can also type `make clean` to clean up any executable or .o files left from previous compilations.
To compile the code documentation execute `make documentation`, this should create a doxigen documentation structure under docs.

<a name="execution"></a>
Running the tool
------

This section describes how to run the system. There is a VERY simple way and a more complicated way.
The simple way is oriented towards casual users of the tool. The more complicated way is part of the legacy code that speech processing people have been using for beamforming audio signals.

<a name="execute_simple"></a>
#### simple (but limited) way

Use the script `do_beamforming.sh` provided in the base directory to apply BeamformIt to all audio files within a certain directory.
The script takes 2 parameters:
- The directory where to find the audio files to be beamformed (it takes all files ending with .wav and .sph)
- The output name "show_id" we want to give to the files
The script will process all files and create a directory with the output files under `./output/$show_id`

<a name="execute_complex"></a>
#### more complicated way

This way allows for the use of a single channels file to process multiple file configurations and to have all audio files for different configurations in the same directory. This is the way that speech researchers have been using BeamformIt until now (with the difference that some config parameters might have now changed)

To explain the way to run it we follow an example based on the RT06s NIST evaluation for meetings (conference room data).

*Prerequisites*: To run the beamforming we need to have the input files in .sph (sphere) or .wav format, containing one or more channels per file (having a set of files with mixed number of channels is fine). Sometimes a preprocesing is performed to the data prior to beamforming. A usual preprocessing step which has given good results at ICSI is Wiener filtering each individual channel.

*Config files*: There are 2 files that need editing before running the beamforming. For this example they are a config file `cfg-files/RT06s_conf.cfg` and a channels file `cfg-files/channels`.
The config file determines the location of all data and the running parameters, as well as the location of the channels file. The channels file contains information of which individual audio files will be combined into a single output file. The config file is the only mandatory parameter to the executable with `-C` .

Let us see in detail the format and possible parameters in each case:

* `cfg-files/channels`: contains the list of input channels. This file is the one that solely determines how many (and which) channels we are using in the beamforming, therefore this file needs to be different for MDM (multiple distant microphones), ADM (all distant microphones) and other combinations of microphones processed. 
For any given desired beamforming output you need to insert a line in the channels file with the format: `show_id list_of_files`
Where `show_id` must be the same alphanumeric string as the one later passed to the system (through command line or config file) through the parameter `show_id`, identifying which show/ID to process. 
The `list_of_files` is a list of *space-separated* audio files to be beamformed together. The actual paths to each of these files do not necessarily need to be defined inside the channels file (it would be very redundant and would create very long lines). BeamformIt will compose the file names by concatenating each `file_name` in `list_of_files` with the config variable `source_dir` as `source_dir/file_name`.
The total number of channels is obtained automatically from this list of files. Be careful with the inclusion of extra spaces in between files or at the end, as these will be treated as extra files and will cause the system to crash.
The files in the channels file lists can either contain acoustic data for a single channel or data for multiple channels (as is the case of AMI circular arrays data). BeamformIt will process either case without any problem. In the case of multiple channels it will create a set of temporal files containing only 1 channel each. These are created in the output directory and can be deleted by hand after execution.

* `cfg-files/RT06s_conf.cfg`: contains most of the configuration parameters passed to the system. The rest of the parameters are passed through command line arguments and usually refer to parameters changing for each show.
For all parameters there normally is a default value as optimized for the ICSI speaker diarization system and documented in [my thesis](www.xavieranguera.com/phdthesis) (note that names might have changed a bit over the years...). A similar explanation of each parameter is also available running the program with --help. The typical parameters in the configuration file are:

    Parameter | default value | Description
    --- | --- | ---
    **scroll_size** (-r) | 250  |   scrolling size used to apply the delays and to output the signal
    **window_size** (-w) | 500  |   cross correlation computation window size. It starts at each scroll point and extends passed its length.
    nbest_amount | 1         |    amount of maximum points for the xcorrelation taken into account as the possible TDOA.
    do_noise_threshold | 1   |    flag wether to apply an automatic noise thresholding. Possible values are:  0-> do not apply it, 1-> apply a percentage threshold (set by noise_percent), 2-> apply an absolute value threshold (set by min_xcorr).
    noise_percent | 10      |    Percentage of frames with lower xcorr taken as noisy (min 0, max 100) , used when do_noise_threshold=1. 
    min_xcorr | 0.1           |         Absolute value (min 0, max 1) used to threshold the frames xcorr to determine them as noisy,  used when do_noise_threshold=2.
    do_optimum_delays | 0      |   Flag that determines whether to apply any postprocessing to the computed delays before running the sum. 
    do_acoustic_modelling | 0  |   Defines which continuity filter to apply to the TDOA values, none(0), viterbi(1) or a simple filter(2)
    trans_weight_nbest | 25    |  (If do_acoustic_modelling==1) Weights applied to the transition probabilities in the Viterbi TDOA decoding phase. This first weight applies to the monochannel Viterbi (first pass).
    trans_weight_multi | 25    |  (If do_acoustic_modelling==1) Same as before, this applying to the multichannel case (pass 2).
    do_avoid_bad_frames | 1 |    flag wether to use the bad frames in the sum process. If set to 1 it filters out those segments (of size scroll_size) which are of a very poor quality.
    do_compute_reference | 2  |   This parameter selects how to define the reference channel for processing. 0 uses the parameter reference_channel as the reference, 1 computes it with an xcorr-based metric and 2 uses an adaptive reference channel along the meeting
    reference_channel | 0    |   If the previous is set to 0, this determines which channel (ordered according to its order in the channels file) is the reference.
    do_use_uem_file | 0       |   flag wether to use a uem file or not(process all the file)
    do_adapt_weights | 1       |  flag wether to use an adaptative weights scheme or fixed weights
    do_write_sph_files | 1     |  flag wether to output the sph files or just run the system to create the auxiliary files
    full_path | 1              |  This flag determines the way that data is read from the channels file, as explained previously.
    uem_file  |              |     Location of the UEM file (in case it is used)
    **source_dir** (-i)   |      |   Location of the source files, this points to the root directory for files in the channels file to be found
    **channels_file** (-c)  |    |   location of the channels file
    **show_id** (-s)  |        |     Show_id of the show/meeting being processed (for example NIST_20030925-1517)
    **result_dir** (-o)  |      |    Directory where to store the results (doesn't need to exist). Subdirectories will be created in it, one for each meeting processed.
    help (-h)     |           |     produces this help message
    config_file (-C)   |            |    config file to be used
    print_features | 1        |   Prints all the feature values after reading all configuration files
    uem_file    |            |      Optional, insert a NIST UEM file to define the regions to process
    delay_variance | 5      |     maximum allowed variance (+- number of lags) for the delays to consider it comes from the same speaker
    ovl_variance | 20       |      maximum allowed variance when calculating the overlap
    do_indiv_channels | 0     |    flag to output individual channels or not. This option enables the output of an individual file for each channel used in the beamforming where each segment has been delayed and its amplitude adjusted using the computed values. If all channels obtained using this switch were summed-up we would obtaine the standard beamforming output.
    output_format | 0	        |     Selects the audio file output format: 0. same format as the input file; 1. sph 16bit; 2. wav 16bit; 3. wav floating point
    delays_margin | 30       |    +- ms around where we search for the peaks of the autocorrelation. This is an important parameter, as it determines how far away can microphones be. A rough estimation can be obtained by delays_margin * 360 = max meters of distance. It is good to set this parameter as tight as possible for your setup, as the system will then run faster and with less probability to select wrong delays.   
    output_second_wav | 0    	|    Determines wether the system writes out a wav file with the 2nd best delays
    print_level | 2        |      Determines the amount of information printed out by the program in STDOUT. Roughly, the levels have been set according to restrictiveness in the following manner: <ul><li>10. the error messages</li><li>5.  the warning messages</li><li>2.  main program messages indicating general progress</li><li>1.  progress messages from the individual functions</li><li>0.  messages of progress at segment level</li></ul>
    do_output_residual | 0   	|		Flag to write out ONLY the residual between two signals. Does not perform any beamforming. It can be used to test how different two signals are from each other.
    do_compute_skew | 0	|		Flag to compute an initial skew/alignment between input signals. This has been introduced in version 3.45 at the same time that internal predefined skew computation has been eliminated. If processing channels which you expect that have some skew, do apply this parameter. An example are the ICSI meetings in the NIST-RT evals.
    skew_margin | 1000	|		+- ms around where we search for the skew between signals 

Parameters in *bold* are necessary for the system to run and without them it will not start. The letters in parenthesis indicate short ways to refer to the same parameter when used as a command line argument, otherwise you need to use --parameter_name.

* `run-files/run_rt06s` (optional): This file contains the commands to execute all of the meetings in the RT06s (conference room) test set. In this case all parameters are taken by default, only the configuration file location and the show_id are defined in the command line.


<a name="output"></a>
Output files
------

After the system runs (and dumps into stdout a lot of stuff...) an output directory is created in result_dir with the same name as `show_id` (in here `NIST_20051024-0930` as an example) and with some/all of the files inside (depending on the config parameters being used):

- `NIST_20051024-0930.sph` : sphere file containing the beamformed output signal, the main output of the system

- `NIST_20051024-0930_2.sph` : a second sphere file (will be empty if output_second_best_sph=0) containing the beamforming output using the 2nd best TDOA delays.

- `NIST_20051024-0930.del` : The best TDOA delays used in computing the beamforming. It contains a delay for each scroll segment and their GCC-PHAT values

- `NIST_20051024-0930.del2` : Same as before, but for the 2nd best delays

- `NIST_20051024-0930.ovl` : Overlap data computed using a xcorrelation-based algorithm where 0/1 are output for each scroll frame. A 1 indicated an overlap region.

- `NIST_20051024-0930.weat` : weights used for each channel at each output segment.

- `NIST_20051024-0930.info` : general information from the system.


Known limitations
------

The BeamformIt software has been tested for up to 128 channels in parallel, but in theory it can hold as many as memory can manage where you run it. To increase the number of channels you just need to change the variable MAXNUMCH in src/global.h

The current version is independent on the number of channels per input file, the framerate and its resolution.
If the input files defined in the "channels" file contain more than 1 channel, an algorithm is executed internally to split the data into individual files, storing them into the output directory in WAV format.
All these files are kept in the directory after the execution, it is the user's choice to delete them if disk space is a constraint.

<a name="cite"></a>
How to Cite
------

If you use this software in your research please reference it as:

"Acoustic beamforming for speaker diarization of meetings", Xavier Anguera, Chuck Wooters and Javier Hernando, IEEE Transactions on Audio, Speech and Language Processing, September 2007, volume 15, number 7, pp.2011-2023.
```@ARTICLE{Anguera-IEEE-07,
  author = {X. Anguera and C. Wooters and J. Hernando},
  title = {Acoustic beamforming for speaker diarization of meetings},
  journal = {IEEE Transactions on Audio, Speech, and Language Processing},
  year = {2007},
  volume = {15},
  pages = {2011-2021},
  number = {7},
  month = {September}
}
```

And also a (very exhaustive) description of the system, with speaker diarization experiments,  
in my [PhD thesis](http://www.xavieranguera.com/phdthesis/):

Xavier Anguera, "Robust Speaker Diarization for Meetings", PhD thesis, Technical University of Catalonia, 2006

```@PHDTHESIS{Anguera06phdthesis,
  author = {Xavier Anguera},
  title = {{PhD Thesis: Robust Speaker Diarization for Meetings}},
  school = {Universitat Politecnica de Catalonia},
  year = {2006}
}
```

<a name="changelog"></a>
Change log
------

This contains the main changes in the tool. This replaces the versioning system we used until version 3.51. For a more detailed changelog refer to the Github commit logs.

- Version 1.0: Initial public version, includes some documentation, a "clean" code and RT06s conf. room example scripts.
- Version 1.1:

    * Source code now has its own directory, making the structure a bit more organized
    * Some warnings and obsolete functions were eliminated in the source code.
    * Documentation now can be compiled directly from the Makefile
    * Other bits and pieces.

- Version 2.0: Some fundamental improvements were made from verion 1.1, therefore I decided to change a major number

    * Solved some memory leaks potentially dangerous.
    * Added a fourth library requirement: libsamplerate and use it to automatically convert to 16k all input files with different samplerate.
    * Made the code independent on the input format data. When it is of different format than expected it automatically converts it.
    * All temporal wavefiles are format WAV in floating point, little endianness.
    * Internal rework of all variables to handle acoustic data in floating point precision, therefore accepting a broader range of input data.
    * Output signal weighting_factor and output normalization changed to accomodate the floating point data.
    * Changed several libraries to link dynamically when possible.
    * Input parameters functions made robust to errors in the parameters passed in. In the previous version the system crashed without an error message if a parrameter was misspelled.
    * Added running and configuration scripts for RT06s lecture room data, changed the names of the conf. room data to accomodate.
    * Variable min_xcorr used some hard-coded values according to each meeting. These have been eliminated.
    * Solved a bug when exctracting multiple channels into independent files. This function has been sped-up by using buffers.
    * Solved a bug that caused the signal to groud forever when the input signal was 0 for one analysis window.
    * Implemented different levels of printout information to accomodate different uses.

- Version 2.1: Some bugfixes solved from 2.0
	* When  outputing the weights into the .weat files it was sometimes outputing 0's when it meant equal weight
	* Also with the final weights, sometimes they were adulterated from their real value. Now they do have them (this should change the output wavefile)
	* When the showID was not following that of NIST meetings, or it was too long, the system encountered serious problems.

- Version 3.0:
    * Another major rework is made, converting all the code into a more C++ style. NO longer a huge functions.cpp file
    * Added a .vcproj project to open Beamformit as a project in Visual Studio.

- Version 3.1:
    * Implemented the variable reference channel selection algorithm, which allows the system to switch reference channel to use that with highest quality according to xcorr metrics. You can enable it with do_compute_reference=2
    * Some parameters tuning is performed using signal quality tools and some default parameters are changed. In the cfg-files two versions are now present for the example meeting sets, with the old and new parameters.

- Version 3.3:
    * Some major bugs are fixed which did not allow to run do_compute_reference=2 with acoustic optimization

- Version 3.4:
	* Solved the bugFixes also included in version 2.1
	* Solved a serious bug in the acoustic modeling of delays (postprocessing step)
	* Solved some Valgrind issues which seemed not dangerous but might had been so in certain cases.
	* Added an option in the menu to compare two output files.
	* Reverted default parameters to those extensively tested in my PhD thesis.
	
- Version 3.4.1:
	* Solved a bug when reading multichannel files regarding allocation of memory, which caused a segmentation fault.
	* expanded the column width of the help printout
	* Changed the option "do_floating_point" to "output_format" to use it as output sound format selection parameter.
	* Implemented standard WAV file sound output. Others can easily added under request 

- Version 3.5 (June 2014):
	* Eliminated the dependency with Boost-command_options library
		A new class has been included (parse_options) to parse command line and config file parameters. Same format as before is used.
	* Eliminated the dependency with libsamplerate library
		In previous releases any input file with sampling rate different than 16KHz were changed to 16KHz before continuing. This does not correspond to any
		need in the tool, as all internal parameters has been adapted to work with the actual sampling rate. We therefore raise this limitation and thus stop using this library.
	* Usage of parameter "source_dir" has been changed. It now contains a directory that will be appended (preceding) any file info obtained from the channels file
	* Parameter "full_path" has been eliminated. All cases now are considered full path
	* Parameter "result_dir" indicates the final directory where all results are stored (before it was combined with "show_id" to create such directory) 
	* Generalization of the cross-correlation maxima search parameters ("delays_margin" parameter).
		In previous versions of BeamformIt the parameters "delays_margin" and "skew_margin" were set inside the program. These determine how far in time the time-delay of arrival is allowed to be.
		These were set for the NIST-RT meetings as: 10ms (LDC, NIST), 15ms (CMU), 20ms (ICSI) and 20ms others.
		Considering a speed of sound at 360m/s this means that 10ms -> 3.6m max distance between microphones.
		In order to make BeamformIt more generic we eliminated these hardcoded values and set it by default in "delays_margin" to 30ms (approx 10m max distance between microphones)
		For meetings data (all microphones are close to each other) smaller values are possible.
		We also created parameter "skew_margin" that look for bigger delays between channels that is set initially and kept as a constant offset for the rest of the processing. This was previously used 
		for ICSI meetings in NIST-RT as these had a bug that started recording in different channels at different times. default value is 1sec.
	* A new parameter (do_compute_skew) is added to select wether to compute the skew/alignment between input signals.
		Previously this was automatically done for ICSI meetings only. If still processing these meetings you will need to turn it on it manually.
	* Renamed "output_second_best_sph" to "output_second_wav" which seems more adequate for the variable
	* The output of the beamforming is renamed to the show name (we eliminate the "_seg" suffix). Same for delay files
	* A script "do_beamforming.sh" has been added to ease the beamforming process for casual users of the tool

- Version 3.51 (July 2014):
  This version contains mostly bugfixes of 3.5. This has been tested to perform equally to previous versions, while incorporating all the advantages of 3.5
  * Added an AMI config and channels file example, for processing the AMI database
  * Reduced the number of completion percentage log output, which was cluttering .log files when piping the output to a log. 

