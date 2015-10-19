/*
	Set of support classes  
*/

#include <stdio.h>
#include "support.h"

/*!
  Default constructor
*/
LogInfo::LogInfo(void)
{
	m_printLevel = 1; //by default
	char tmp_string[1024];sprintf(tmp_string, "No messages to report\n"); //initialized to avoid warnings
}

/*!
  Default destructor
*/
LogInfo::~LogInfo()
{
}

/*! 
  Print manager to output data in STDOUT, uses the level determined by configuration to decide wether to print

  \param pstring is the data to be printed
  \param level is the level considered wether to print the information
  \return 0 if did not print anything, 1 if printed
*/

int LogInfo::print_this(string pstring, int level)
{
  if(level >= m_printLevel)
  {
    printf("%s",pstring.c_str());
    fflush(stdout);
    return(1);
  }

  return(0);

}

/*!
  Wrapper for the print manager to pass it a char* instead
  \param pstring is the data to be printed
  \param level is the level considered wether to print the information
  \return 0 if did not print anything, 1 if printed
*/
int LogInfo::print_this(const char* pstring, int level)
{
  string tmpString(pstring);
  return(print_this(tmpString, level));
}

/*!
  Debug function to print the parameters of choice taken from the config structure.
  \param config Is the structure that holds the system configuration parameters
  \return Nothing
 */
void LogInfo::print_features(Configuration &config)
{
  //prints all the features that can be selected via the configuration input, and their associated value
  //print all the parameters:
  print_this("Parameters settings(some of them):\n Mandatory parameters:\n", 1);
	char tmp_string[1024];  
	sprintf(tmp_string,"  scroll_size: %d\n", config.rate); print_this(tmp_string, 1);
  sprintf(tmp_string,"  window_size: %d\n", config.window); print_this(tmp_string, 1);
  sprintf(tmp_string,"  source_dir: %s\n", config.AUDIOROOT.c_str()); print_this(tmp_string, 1);
  sprintf(tmp_string,"  show_id: %s\n", config.SHOWNAME.c_str()); print_this(tmp_string, 1);
  sprintf(tmp_string,"  result_dir: %s\n", config.RESULTPATH.c_str()); print_this(tmp_string, 1);
  sprintf(tmp_string,"  channels_file: %s\n", config.CHANNELSFILE.c_str()); print_this(tmp_string, 1);

  print_this(" Optional parameters:\n", 1);
  sprintf(tmp_string,"  uem_file: %s\n", config.UEMFILE.c_str());  print_this(tmp_string, 1);
  sprintf(tmp_string,"  delay_variance: %d\n", config.delay_variance); print_this(tmp_string, 1);
  sprintf(tmp_string,"  ovl_variance: %d\n", config.ovl_variance); print_this(tmp_string, 1);
  sprintf(tmp_string,"  do_optimum delays: %d\n", config.DO_OPTIMUM_DELAYS); print_this(tmp_string, 1);
  sprintf(tmp_string,"  do_indiv_channels: %d\n", config.INDIV_CHANNELS); print_this(tmp_string, 1);
  sprintf(tmp_string,"  output_format: %d\n", config.OUT_FORMAT); print_this(tmp_string, 1);
  sprintf(tmp_string,"  do_xcorr_reference: %d\n", config.COMPUTE_REFERENCE); print_this(tmp_string, 1);
  sprintf(tmp_string,"  nbest_amount: %d\n", config.nbest_amount); print_this(tmp_string, 1);
  sprintf(tmp_string,"  reference_channel: %d\n", config.reference_channel); print_this(tmp_string, 1);
  sprintf(tmp_string,"  Transition nbest probability: %f\n", config.trans_weight_nbest); print_this(tmp_string, 1);
  sprintf(tmp_string,"  Percentage of noisy xcorr values: %f\n", config.noise_percent); print_this(tmp_string, 1);
}
