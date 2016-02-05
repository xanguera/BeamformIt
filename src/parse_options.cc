/*
  parse_options - Xavier Anguera, 2014
  Class used to parse the input parameters and fill the config class instance with command line and input 
  Parameters can be entered via command line or via config file. When the same parameter is entered in both ways, the command ine option prevails.
  Parameters in the config file have to follow the format "name = value". Anything after '#' is understood as a comment.
  Command line parameters can be input as "--name/-n value" or "--name/-n=value" (where "--name/n" means that some parameters have a short and long version available)  
*/

#include "parse_options.h"

//parses the input parameters
void parse_options::parse_params(int argc, char *argv[])
{
    //register all parameters
    //mandatory parameters:
    add_param("scroll_size", "r", "int", &((*m_config).rate), true, "250", "scroll size (ms for which the same delay is applied)");
    add_param("window_size", "w", "int", &((*m_config).window), true, "500", "Computation window (window in ms where each delay is computed)");
    add_param("show_id", "s", "string", &((*m_config).SHOWNAME), true, "", "The show ID that we want to compute");
    add_param("result_dir", "o", "string", &((*m_config).RESULTPATH), true, "", "Output files base directory");
    add_param("channels_file", "c", "string", &((*m_config).CHANNELSFILE), true, "", "Channels file defining what individual files to beamform");

    //optional parameters
    add_param("source_dir", "i", "string", &((*m_config).AUDIOROOT), false, "", "If set, it gets appended to audio files obtained from th channels file");
    add_param("help", "h", "flag", NULL, false, "false", "produces this help message");
    add_param("config_file", "C", "string", &((*m_config).m_config_file), false, "", "config file to be used");
    //add_param("full_path", "", "int", &((*m_config).full_path), false, "1", "Switch to select how the channels file is read: 0)");
    add_param("print_features", "", "int", &((*m_config).print_features), false, "1", "Prints all the feature values after reading all configuration files");
    add_param("uem_file", "", "string", &((*m_config).UEMFILE), false, "", "Optional, insert a NIST UEM file to define the regions to process");
    add_param("delay_variance", "", "int", &((*m_config).delay_variance), false, "5", "maximum allowed variance (+- number of lags) for the delays to consider it the same speaker");
    add_param("ovl_variance", "", "int", &((*m_config).ovl_variance), false, "20", "maximum allowed variance when calculating the overlap");
    add_param("do_optimum_delays", "", "int", &((*m_config).DO_OPTIMUM_DELAYS), false, "1", "flag whether to apply delays postprocessing[0]");
    add_param("do_acoustic_modelling", "", "int", &((*m_config).DO_ACOUSTIC_MODELLING), false, "1", "When doing delays postprocessing, defines which N-best selection to apply: 0->none, 1->viterbi or a 2->simple filter[0]");
    add_param("do_indiv_channels", "", "int", &((*m_config).INDIV_CHANNELS), false, "0", "flag to output individual channels or not");
    add_param("output_format", "", "int", &((*m_config).OUT_FORMAT), false, "0", "Selects the audio file output format: 0. same format as the input file; 1. sph 16bit; 2. wav 16bit; 3. wav floating point");
    add_param("do_compute_reference", "", "int", &((*m_config).COMPUTE_REFERENCE), false, "1", "Sets the method to select the reference channel: 0-> set by reference_channel config variable, 1-> set at start by xcorr, 2-> dynamically set by xcorr");
    add_param("reference_channel", "", "int", &((*m_config).reference_channel), false, "0", "Sets the reference channel number, overwrites any other computation");
    add_param("min_xcorr", "", "float", &((*m_config).min_xcorr), false, "0.1", "minimum absolute value of the GCC-PHAT xcorrelation to consider a delay when doing fixed noise thresholding");
    add_param("do_noise_threshold", "", "int", &((*m_config).DO_NOISE_THRESHOLD), false, "", "Flag wether to apply an automatic noise thresholding or not, and either if it is fixed of a percent");
    add_param("noise_percent", "", "float", &((*m_config).noise_percent), false, "10", "Percentage of frames with lower xcorr taken as noisy");
    add_param("delays_margin", "", "int", &((*m_config).margin), false, "30", "+- ms around where we search for the peaks of the autocorrelation (30ms ~= 11 meters max between microphones)");
    add_param("skew_margin", "", "int", &((*m_config).skewMargin), false, "1000", "+- ms around where we search for the skew between signals");
    add_param("nbest_amount", "", "unsigned int", &((*m_config).nbest_amount), false, "4", "Number of best xcorr values to consider[1]");
    add_param("do_avoid_bad_frames", "", "int", &((*m_config).DO_AVOID_BAD_FRAMES), false, "0", "Flag to use the bad frames (according to xcorr) during the sum process");
    add_param("do_use_uem_file", "", "int", &((*m_config).DO_USE_UEM_FILE), false, "1", "Flag to use a uem file or process the entire file");
    add_param("do_adapt_weights", "", "int", &((*m_config).DO_ADAPT_WEIGHTS), false, "1", "Flag to adapt the weights according to xcorr or fix the weights");
    add_param("do_write_sph_files", "", "int", &((*m_config).DO_WRITE_SPH_FILES), false, "1", "Flag to write out the sph files or just the aux files");
    add_param("trans_weight_nbest", "", "float", &((*m_config).trans_weight_nbest), false, "25", "transition probs weight");
    add_param("trans_weight_multi", "", "float", &((*m_config).trans_weight_multi), false, "25", "transition probs weight");
    add_param("output_second_wav", "", "int", &((*m_config).output_second_best_sph), false, "0", "Determines wether the system writes out a second audio file obtained from the 2nd best delays");
    add_param("print_level", "", "int", &((*m_config).m_printLevel), false, "1", "Strictness level applied when printing informations. 0 prints everything, inf. prints nothing");
    add_param("do_output_residual", "", "int", &((*m_config).DO_OUTPUT_RESIDUAL), false, "0", "Flag to write out ONLY the residual between two signals. Does not perform any beamforming");
    add_param("do_compute_skew", "", "int", &((*m_config).DO_COMPUTE_SKEW), false, "0", "Flag to compute an initial skew/alignment between input signals");
    add_param("do_output_overlap", "", "int", &((*m_config).DO_OUTPUT_OVERLAP), false, "0", "Flag indicating whether to compute overlaps from the TDOA's (VERY experimental feature)");
    add_param("do_signal_quality", "", "int", &((*m_config).DO_SIGNAL_QUALITY), false, "0", "flag whether to compute signal quality from the resulting audio signal");
    add_param("bad_frames_ratio", "", "int", &((*m_config).badFramesRatio), false, "10", "Sets the threshold to separate bad from good frames, the higher the value, the less bad frames");
    add_param("ref_adapt_ratio", "", "int", &((*m_config).refAdaptRatio), false, "0.05", "Dynamic reference channel selection adaptation ratio");
    add_param("do_delays_padding", "", "int", &((*m_config).DO_DELAYS_PADDING), false, "1", "Flag to add extra delays at the end of the computed ones to try filling all the acoustic data, which were not computed due to a large analysis data");
    add_param("extra_delays_padding", "", "int", &((*m_config).extra_delays_padding), false, "10", "Number of extra padding delays appended in addition to the regular padding to the end of the computed delays. Relevant to ensure the number of delays exceeds the length of the file in some applications");

  //parse all input params from the command line
  if(processCommandLine(argc, argv) == false)
  {
      exit(-1); //TODO: actually return an error
  }

  //process, if available, the config file
  string configFile;
  if(m_params["config_file"]->m_externallySet)
  {
      configFile = *(string*)(m_params["config_file"]->m_value);
      printf("Processing config file %s\n", configFile.c_str());
      if(!configFile.empty())
          processConfigFile(configFile);
  }

  //check whether the user requested help
  if(m_params["help"]->m_externallySet)
  {
      print_help();
      exit(0);
  }

  //check whether all mandatory parameters have been filled
  for(map<std::string, Params*>::iterator iter = m_params.begin(); iter != m_params.end(); iter++)
  {
      //printf("Checking %s\n", iter->second->m_longID.c_str());
      if(iter->second->m_mandatory && !iter->second->m_externallySet)
      {
          printf("ERROR: parameter %s has not been set and it is mandatory. Run with -h for help\n", iter->second->m_longID.c_str());
          exit(-1);
      }
  }
  
}


//process the command line options
bool parse_options::processCommandLine(int argc, char *argv[])
{
    int count=1;
    while(count < argc)
    {
      //get the parameter name
      string paramID(argv[count]);
      //The parameter name can either be in short form "-X" or in long form "--XXX"
      //eliminate the "-" and "--" from the ID
      string::iterator end = paramID.begin()+1;
      if(paramID[1] == '-')
        end++;
      paramID.erase(paramID.begin(), end);

      //I accept two forms of input parameters: "--name=value" or "--name value"
      //Look for the existance of "=" in the paramID
      string paramValue;
      bool equalType = false; //true if I find a '='
      std::size_t found=paramID.find('=');
      if (found!=std::string::npos)
      {
        //we have in this case the format name=value
        paramValue = paramID.substr(found+1);
        paramID.erase(found);
        equalType = true;
      }      
      
      //find the parameter long version and whether it exists
      if(m_options.find(paramID) == m_options.end())
      {
          printf("ERROR: input option %s is not available\n", paramID.c_str());
          print_help();
          exit(0);
      }
      string paramName = m_options[paramID]; //get the right key

      //check whether it is a flag or a key-value pair
      if(m_params[paramName]->m_type == "flag")
      {
          paramValue="true";
      }
      else
      {
        if(!equalType)
        {
          //get the parameter value
          if(count+1 >= argc)
          {
              printf("ERROR: option %s does not have a value\n", paramID.c_str());
          }
          count++;
          paramValue.assign(argv[count]);
        }
      }

      //test
      //printf("Registering <%s> -> <%s> with value <%s>\n", paramID.c_str(), paramName.c_str(), paramValue.c_str());

      //store the parameter value
      if(m_params.find(paramName) != m_params.end())
      {
          //the parameter exists
          if(m_params[paramName]->m_externallySet)
              printf("WARNING: parameter <%s> with value <%s> has already been defined, prior one will be overwriten\n", paramID.c_str(), paramValue.c_str());
          m_params[paramName]->set_value(paramValue);
      }
      else
      {
          //the parameter does not exist
          printf("ERROR: input option %s is not available\n", paramID.c_str());
          print_help();
          exit(0);
      }
      count++;
    }
    return true;
}


//process config file
void parse_options::processConfigFile(string configFile)
{
    FILE * pFile;
    pFile = fopen (configFile.c_str() , "r");
    if (pFile == NULL)
    {
      printf("ERROR opening config file <%s>", configFile.c_str());
      exit(-1);
    }
    else
    {
      char buffer[512];
      while ( ! feof (pFile) )
      {
        if ( fgets(buffer , 512 , pFile) == NULL ) break;

        //convert line feed by \0
        buffer[strlen(buffer)-1] = '\0';
        
        //process the input line
        if(buffer[0] == '#' || buffer[0] == '\0') continue; //we have a comment or empty line
        
        //we look for an option of the format "optionID = paramValue"
        int i=0;
        while(buffer[i] != ' ' && buffer[i] != '=') //look for the end of the optionID
          i++;
        string paramID(buffer, i);

        while(buffer[i] == ' ' || buffer[i] == '=') //eliminate the middle characters
          i++;

        int startValue=i;
        while(buffer[i] != ' ' && buffer[i] != '#' && buffer[i] != '\0') //find the end of the second
          i++;
        string paramValue(&(buffer[startValue]), (i-startValue));

        //register that ID-value pair
        //store the parameter value
        string paramName = m_options[paramID];

        //test
        //printf("Registering <%s> -> <%s> with value <%s>\n", paramID.c_str(), paramName.c_str(), paramValue.c_str());

        if(m_params.find(paramName) != m_params.end())
        {
            //the parameter exists. I give priority to the command line
            if(m_params[paramName]->m_externallySet)
                printf("INFO: parameter <%s> with value <%s> has already been defined with value <%s>. Prior value is kept\n", paramID.c_str(), paramValue.c_str(), paramValue.c_str());
            else
                m_params[paramName]->set_value(paramValue);
        }
        else
        {
            //the parameter does not exist
            printf("ERROR: input option %s is not available\n", paramID.c_str());
            print_help();
            exit(0);
        }
      }
      fclose (pFile);
    }
}


//adds a new parameter into the system
int parse_options::add_param(string longID, string shortID, string type, void* destination, bool mandatory, string defaultValue, string param_desc)
{
  Params* tmpParam = new Params;
  tmpParam->m_longID = longID;
  if(shortID != "")
      tmpParam->m_shortID = shortID;
  tmpParam->m_value = destination;
  tmpParam->m_mandatory = mandatory;
  tmpParam->m_default = defaultValue;
  tmpParam->m_param_desc = param_desc;
  tmpParam->m_type = type;
  //we set the value
  tmpParam->set_value(defaultValue);
  tmpParam->m_externallySet = false; //initially not set


  //we store the value into the long map, and keep a mapping in case short params are used
  //check whether the same ID's have been used before
  if(m_params.find(longID) != m_params.end())
  {
      printf("ERROR: An input parameter with the same ID as %s has been registered before\n", longID.c_str());
      exit(-1);
  }
  m_params[longID] = tmpParam;
  //m_params[longID] = tmpParam;
  m_options[longID] = longID;
  if(shortID != "")
      m_options[shortID] = longID;

  return(true);
}

void parse_options::print_help()
{
    //iterate through the map containing defined options
    printf("BeamformIt options:\n");
    printf("(Options with * are mandatory)\n");
    for(map<std::string, Params*>::iterator iter = m_params.begin(); iter != m_params.end(); iter++)
    {
        printf("%s, %s [ = %s ]",iter->second->m_longID.c_str(), iter->second->m_shortID.c_str(),iter->second->m_default.c_str());
        if(iter->second->m_mandatory) printf(" * ");
        printf("\t\t%s\n", iter->second->m_param_desc.c_str());
    }
}



/*
 *Inserts an input parameter into the destination, depending on the type
 */
bool Params::set_value(string value)
{

    if (m_type == "string")
    {
        *(std::string*) m_value = value;
    }
    else if (m_type == "int" || m_type == "unsigned int")
    {
        *(int*) m_value = atoi(value.c_str());
    }
    else if (m_type == "unsigned long")
    {
        *(unsigned long*) m_value = atoi(value.c_str());
    }
    else if (m_type == "float")
    {
        *(float*) m_value = atof(value.c_str());
    }
    else if (m_type == "double")
    {
        *(double*) m_value = atof(value.c_str());
    }
    else if (m_type == "bool" || m_type == "flag")
    {
        //we use the bool/flag for the help, in which case I will have NULL
        if(m_value != NULL)
        {
            if(value == "true")
                *(bool*) m_value = true;
            else if(value == "false")
                *(bool*) m_value = false;
            else
            {
                printf("ERROR: bool type can only accept true/false\n");
                return false;
            }
        }
    }
    else
    {
        printf("ERROR: type %s is not supported\n", m_type.c_str());
        return false;
    }

    //we tell that the parameter has been set
    m_externallySet = true;

    return true;
}

