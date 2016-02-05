#include <stdlib.h>
#include <cstring>
#include "fileinout.h"
//#include "samplerate.h"

//char tmp_string[1024];

/*!
  Default desctructor
*/
FileInOut::~FileInOut()
{
  //cleans all allocated memory
  for(int i=0; i<m_numFiles; i++)
  {
    delete[] m_channels[i];
  }
	delete[] m_fileOut;
}

/*!
  Constructor setting the config param
  \param config pointer to the config class
*/
FileInOut::FileInOut(Configuration* config)
{
	m_config = config;
}

/*
	Functions for input-output of data
*/

/*!
  Initialization of the class
  \return nothing
*/
void FileInOut::init()
{
}

/*!
  It first gets the information about the files being processed from the channels file. It reads in the
  names of the input files, does not extract the channels information yet

  It then opens input files to be processed and obtains the number of channels to be processed, allowing for each file to
  contain multiple channels each. It such case, it extracts all channels in the output directory as individual channel files,
  and opens those.

  \return Number of channels read

*/

int FileInOut::Open_Input_Channels()
{

    ///// First looks into the channels file and extracts the file list /////

    m_channels.clear();
    FILE *chanfd;
    chanfd = fopen((*m_config).CHANNELSFILE.c_str(), "r");
    if(chanfd == NULL)
    {
        char tmp_string[1024]; sprintf(tmp_string, "Error Opening file %s\n",(*m_config).CHANNELSFILE.c_str());
        m_log.print_this(tmp_string, 10);
        exit(1);
    }

    int count=0;
    int scroll=0;
    int num_files = 0; //counts the number of files in the channels file

    char tmp_string[1024]; sprintf(tmp_string, "Extracting channel information for show %s\n", (*m_config).SHOWNAME.c_str());
    m_log.print_this(tmp_string, 1);

    char line[4096];
    while(fgets(line, 4096, chanfd) != NULL)
    {
        //check if it's the line we want
        if(!strncmp(line, (*m_config).SHOWNAME.c_str(), (*m_config).SHOWNAME.length()))
        {
            count = (*m_config).SHOWNAME.length(); //we start after the show ID
            while(line[count] != '\0' && line[count] != '\n' && line[count] != '\r') //'\r' for windows files
            {
                count++; //from the space to the string
                //new file ending
                char *tmp_channelName = new char[2048];
                m_channels.push_back(tmp_channelName);
                scroll = count; //the beginning of this string
                while(line[count] != ' ' && line[count] != '\0' && line[count] != '\n' && line[count] != '\r')
                {
                    m_channels[num_files][count-scroll] = line[count];
                    count++;
                }
                m_channels[num_files][count-scroll] = '\0';

                char tmp_string[1024]; sprintf(tmp_string, "File %d: %s\n",num_files, m_channels[num_files]);
                m_log.print_this(tmp_string, 1);
                num_files++;
            }
        }
    }
    fclose(chanfd);
    //check whether we did find the channel info
    if(num_files == 0)
    {
      printf("[ERROR]: No information on showID %s was found in channels file %s\n", (*m_config).SHOWNAME.c_str(), (*m_config).CHANNELSFILE.c_str());
      exit(-1);
    }

    m_numFiles = num_files;

    ////// Then opens each file and gets its channel#, if it is >1, it creates temporal files
    ////// and stores the individual data

    char file_name[1024];
    int num_channels = 0;
    for(int count=0; count< m_numFiles; count++)
    {
        //if((*m_config).full_path)
        //{
        sprintf(file_name,"%s/%s", (*m_config).AUDIOROOT.c_str(), m_channels[count]);
        printf("Filename: %s\n", file_name);
        //}
        //else
        //{
        //  sprintf(file_name,"%s%s/%s%s.sph",(*m_config).AUDIOROOT.c_str(), (*m_config).SHOWNAME.c_str(), (*m_config).SHOWNAME.c_str(), m_channels[count]);
        //}

        //retrieve the extension of the input file (we expect it to have length 3)
        strncpy(input_file_extension, &file_name[strlen(file_name)-3], 3);
        input_file_extension[3]='\0';

        //printf("File name <%s> Extension is <%s>\n", file_name, input_file_extension);

        file_info.format = 0; //mandatory when opening for reading
        inputfd[num_channels] = sf_open(file_name, SFM_READ, &(file_info));

        if(inputfd[num_channels] == NULL)
        {
        /*
        //sometimes we have a sph file with .wav extension. Here is a dirty way to check if that's the problem
      if(!(*m_config).full_path)
      {
        //we construct the file_name2 as the initial, substituting sph by wav
        char file_name2[1024];
        strcpy(file_name2, file_name);
        strcpy(&file_name2[strlen(file_name2)-3], "wav");

        file_info.format = 0; //mandatory when opening for reading
        inputfd[num_channels] = sf_open(file_name2, SFM_READ, &(file_info));

        if(inputfd[num_channels] == NULL)
        {
          char tmp_string[1024]; sprintf(tmp_string, "ERROR: Input file %s open error\nERROR: Input file %s open error\n", file_name, file_name2);
          m_log.print_this(tmp_string, 10);
          exit(1);
        }
        else
        {
          //opened correctly, replace the correct file name
          strcpy(file_name, file_name2);
        }
      }
      else
      {
      */
        char tmp_string[1024];
        snprintf(tmp_string,
          sizeof(tmp_string),
          "ERROR: Input file %s open error: %s\n",
          file_name,
          sf_strerror(NULL));
        m_log.print_this(tmp_string, 10);
        exit(1);
        //}
        }
        else
        {

            char tmp_string[1024]; sprintf(tmp_string, "File %s (%d channels) open\n", file_name, file_info.channels);
            m_log.print_this(tmp_string, 1);
        }

        //check that the sampling rate is the same for all files. I do not support different sampling rates
        if(count == 0)
            m_sampleRate = file_info.samplerate;
        else
            if(m_sampleRate != file_info.samplerate)
            {
                printf("[ERROR] file %s has sampling rate %d, different from previous files with sampling rate %ld.\n"
                       "BeamformIt can only run when all files have the same sampling rate.\n"
                       "You can fix it using sox\n", file_name, file_info.samplerate, m_sampleRate);
                exit(1);
            }

        //save the original input format
        file_input_info = file_info;

        /*
    //check the sampling rate: if modified, converts to sph, 16k
    //TODO: need to make BeamformIt independent of the sampling rate!!!!
    if(file_info.samplerate != 16000)
    {
        char tmp_string[1024]; sprintf(tmp_string, "Samplerate is %d, Beamformit works optimally for 16K files\n",file_info.samplerate);
        m_log.print_this(tmp_string, 1);
        //create the name
        char file_out[512];
        sprintf(file_out,"%s%s_f%d_16k.wav",(*m_config).RESULTPATH.c_str(), (*m_config).SHOWNAME.c_str(), count);
        sprintf(tmp_string, "Creating wav file (float point) at 16k (little endian) in %s\n", file_out);
        m_log.print_this(tmp_string, 1);

        int number_samples = file_info.frames;

        int samplerate_in = file_info.samplerate;


        SNDFILE* outputfd;
        file_info.samplerate = 16000; //force to be 16K
        int samplerate_out = file_info.samplerate;

        if(!sf_format_check(&(file_info)))
        {
            m_log.print_this("WARNING: Something is wrong with the file description\n", 5);
        }

        //change format to WAV float from whatever was before
        file_info.format &= ~SF_FORMAT_TYPEMASK; //put to 0 the format bits
        file_info.format |= SF_FORMAT_WAV; //setting them to NIST SPH
        file_info.format &= ~SF_FORMAT_ENDMASK; //set endianness
        file_info.format |= SF_ENDIAN_CPU;
        file_info.format &= ~SF_FORMAT_SUBMASK;
        file_info.format |= SF_FORMAT_FLOAT;

        outputfd = sf_open(file_out, SFM_WRITE, &(file_info));
        if(outputfd == NULL)
        {
            char tmp_string[1024]; sprintf(tmp_string, "ERROR: output file %s open error\n", file_out);
            m_log.print_this(tmp_string, 10);
            exit(1);
        }

        float input_data[samplerate_in * file_info.channels]; //1s worth of data
        float output_data[samplerate_out * file_info.channels]; //1s worth of data

        //create a libsamplerate converter instance
        SRC_STATE* convert;
        int convert_error;
        SRC_DATA conv_data;
        convert = src_new(SRC_SINC_BEST_QUALITY, file_info.channels, &convert_error) ;
        if(convert == NULL)
        {
            char tmp_string[1024]; sprintf(tmp_string, "Samplerate conversor error: %s\n", src_strerror(convert_error));
            m_log.print_this(tmp_string, 10);
            exit(1);
        }

        conv_data.data_in = input_data;
        conv_data.data_out = output_data;
        conv_data.src_ratio = (double)(samplerate_out)/samplerate_in;

        int read_pointer = 0;
        int write_pointer = 0;
        int end_loop = 0;
        while(read_pointer < number_samples)
        {
            //read data in
            sf_seek(inputfd[num_channels], read_pointer, SEEK_SET);
            int frames_in = sf_readf_float(inputfd[num_channels], conv_data.data_in, samplerate_in);

            if(frames_in < samplerate_in)
            {
                end_loop = 1;
            }

            //convert samplerate
            conv_data.input_frames = (long)(frames_in);
            conv_data.output_frames = (long)(samplerate_out);
            conv_data.end_of_input = end_loop;

            convert_error = src_process(convert, &conv_data);
            if(convert_error != 0)
            {
                char tmp_string[1024]; sprintf(tmp_string, "Conversion error occurred: %s\n", src_strerror(convert_error));
                m_log.print_this(tmp_string, 10);
                exit(1);
            }
            read_pointer += conv_data.input_frames_used;

            //write the out data
            sf_seek(outputfd, write_pointer, SEEK_SET);
            sf_writef_float(outputfd, conv_data.data_out, conv_data.output_frames_gen);
            write_pointer += conv_data.output_frames_gen;

        }

        src_delete(convert) ;
        sf_close(outputfd);

        //finally, close the current input file and reopen the new one
        sf_close(inputfd[num_channels]);
        file_info.format = 0; //mandatory when opening for reading
        inputfd[num_channels] = sf_open(file_out, SFM_READ, &(file_info));

    }
    */


    ////////////////

    //check the number of channels
    if(file_info.channels > 1)
    {
      //too many channels, we need to create individual files
      char tmp_string[1024]; sprintf(tmp_string, "Number of channels for file %s is %d, creating temporal mono-channel files (wav, floating point)\n", file_name, file_info.channels);
      m_log.print_this(tmp_string, 1);
      int number_files = file_info.channels;
      int number_frames = file_info.frames;

      file_info.channels = 1; //create single channel files the rest is the same
      file_info.frames = 0;
      //change format to WAV float from whatever was before
      file_info.format &= ~SF_FORMAT_TYPEMASK; //put to 0 the format bits
      file_info.format |= SF_FORMAT_WAV; //setting them to NIST SPH
      file_info.format &= ~SF_FORMAT_ENDMASK; //set endianness
      file_info.format |= SF_ENDIAN_CPU;
      file_info.format &= ~SF_FORMAT_SUBMASK;
      file_info.format |= SF_FORMAT_FLOAT;

      char file_out[512];
      SNDFILE* tmpfd[number_files];

      for(int i=0; i< number_files; i++)
      {
        //create the name
        sprintf(file_out,"%s%s_f%d_ch%d.wav",(*m_config).RESULTPATH.c_str(), (*m_config).SHOWNAME.c_str(), count, i);
        //open
        char tmp_string[1024]; sprintf(tmp_string, "Opening %s file to write\n", file_out);
        m_log.print_this(tmp_string, 1);
        tmpfd[i] = sf_open(file_out, SFM_WRITE, &(file_info));
        if(tmpfd[i] == NULL)
        {
          char tmp_string[1024]; sprintf(tmp_string, "output file %s open error\n", file_out);
          m_log.print_this(tmp_string, 10);
          exit(1);
        }
      }

      //write all data(interlaced) to them (use blocks of 10seconds to agilize the passing)
      int buffer_size = 160000;
      float *in_buffer = new float[number_files * buffer_size];
      float **ptr_out_buffer = new float* [number_files];
      float *out_buffer = new float[number_files * buffer_size];

      for(int i=0; i<number_files; i++)
      {
        ptr_out_buffer[i] = out_buffer + buffer_size*i;
      }
      int frames_read;
      int frame_count = 0;
      while(frame_count < number_frames)
      {
        frames_read = sf_readf_float(inputfd[num_channels], in_buffer, buffer_size);
        //deinterleave the data
        for(int i=0; i<frames_read; i++)
        {
          for(int j=0; j<number_files; j++)
          {
            ptr_out_buffer[j][i] = in_buffer[i*number_files+j];
          }
        }
        //write to the output files
        for(int i=0; i<number_files; i++)
        {
          sf_write_float(tmpfd[i], ptr_out_buffer[i], frames_read);
        }
        frame_count += frames_read;
      }

      //clean up
      delete[] in_buffer;
      delete[] out_buffer;
      delete[] ptr_out_buffer;

      //close input and output files
      for(int i=0; i< number_files; i++)
      {
        sf_close(tmpfd[i]);
      }
      sf_close(inputfd[num_channels]);

      //finally, open these files for the program to use
      char file_in[512];
      for(int i=0; i< number_files; i++)
      {
        //create the name
        sprintf(file_in,"%s%s_f%d_ch%d.wav",(*m_config).RESULTPATH.c_str(), (*m_config).SHOWNAME.c_str(), count, i);
        //open
        char tmp_string[1024]; sprintf(tmp_string, "Opening %s file to read\n", file_in);
        m_log.print_this(tmp_string, 1);
        inputfd[num_channels] = sf_open(file_in, SFM_READ, &(file_info));
        if(inputfd[num_channels] == NULL)
        {
          char tmp_string[1024]; sprintf(tmp_string, "Input file %s open error\n", file_in);
          m_log.print_this(tmp_string, 10);
          exit(1);
        }
        num_channels++;
      }
    }
    else
    {
      //only one channel
      num_channels++;
    }

  }

  //we save the file info
  m_frames = file_info.frames;
  m_framesInit = m_frames;
  m_sampleRate = file_info.samplerate;
  m_sampleRateInMs = m_sampleRate/1000;
  m_numCh = num_channels;

  //amount of frames out each iteration
  m_framesOut = (int)((*m_config).rate*m_sampleRateInMs);
  return(num_channels);
}

/*!
  Opens the ouput channels for signal and data output
  \return Nothing
*/

void FileInOut::Open_Output_Channels()
{
    m_fileOut = new char[512];
    char file_out2[512];
    char tmp_string[1024];

    switch((*m_config).OUT_FORMAT)
    {
    case 0:
        //use the same format as the input
        file_info.format = file_input_info.format;
        sprintf(m_fileOut,"%s%s.%s",(*m_config).RESULTPATH.c_str(), (*m_config).SHOWNAME.c_str(), input_file_extension);
        sprintf(file_out2,"%s%s_2.%s",(*m_config).RESULTPATH.c_str(), (*m_config).SHOWNAME.c_str(), input_file_extension);
        break;
    case 1:
        //sph, 16bit
        file_info.format = SF_FORMAT_NIST | SF_FORMAT_PCM_16;
        sprintf(m_fileOut,"%s%s.sph",(*m_config).RESULTPATH.c_str(), (*m_config).SHOWNAME.c_str());
        sprintf(file_out2,"%s%s_2.sph",(*m_config).RESULTPATH.c_str(), (*m_config).SHOWNAME.c_str());
        break;
    case 2:
        //wav file, 16bit
        file_info.format = SF_FORMAT_WAV | SF_FORMAT_PCM_16;
        sprintf(m_fileOut,"%s%s.wav",(*m_config).RESULTPATH.c_str(), (*m_config).SHOWNAME.c_str());
        sprintf(file_out2,"%s%s_2.wav",(*m_config).RESULTPATH.c_str(), (*m_config).SHOWNAME.c_str());
        break;
    case 3:
        //wav, floating point
        file_info.format = SF_FORMAT_RAW | SF_FORMAT_FLOAT;
        sprintf(m_fileOut,"%s%s.wav",(*m_config).RESULTPATH.c_str(), (*m_config).SHOWNAME.c_str());
        sprintf(file_out2,"%s%s_2.wav",(*m_config).RESULTPATH.c_str(), (*m_config).SHOWNAME.c_str());
        break;
    default:
        sprintf(tmp_string, "ERROR: output format %d is not defined\n", (*m_config).OUT_FORMAT);
        m_log.print_this(tmp_string, 10);
        exit(1);
        break;
    }

    //sprintf(tmp_string, "Current settings are %d %d %d\n", file_info.samplerate, file_info.channels, file_info.format);
    //m_log.print_this(tmp_string, 1);

    //open file for 1st best
    sprintf(tmp_string, "Opening %s file\n", m_fileOut);
    m_log.print_this(tmp_string, 1);
    outfd[0] = sf_open(m_fileOut, SFM_WRITE, &(file_info));
    if(outfd[0] == NULL)
    {
        char tmp_string[1024]; sprintf(tmp_string, "ERROR: output file %s open error\n", m_fileOut);
        m_log.print_this(tmp_string, 10);
        exit(1);
    }

    //open file for 2nd best
    if((*m_config).output_second_best_sph)
    {
        sprintf(tmp_string, "Opening %s file\n", file_out2);
        m_log.print_this(tmp_string, 1);
        outfd[1] = sf_open(file_out2, SFM_WRITE, &(file_info));
        if(outfd[1] == NULL)
        {
            char tmp_string[1024]; sprintf(tmp_string, "ERROR: output file %s open error\n", file_out2);
            m_log.print_this(tmp_string, 10);
            exit(1);
        }
    }

    //we open the individual output channels
    if((*m_config).INDIV_CHANNELS)
    {
        for(int count=0; count< m_numCh; count++)
        {
            //get the basename from the channel files
            string tmpFile(m_channels[count]);
            unsigned found = tmpFile.rfind("/");
            if (found!=std::string::npos)
                tmpFile.erase(tmpFile.begin(), tmpFile.begin()+found);

            sprintf(m_fileOut,"%s/%s",(*m_config).RESULTPATH.c_str(), tmpFile.c_str());
            printf("We will output to file %s\n", m_fileOut);
            //sprintf(m_fileOut,"%s%s%s.sph",(*m_config).RESULTPATH.c_str(), (*m_config).SHOWNAME.c_str(), m_channels[count]);

            outputfd[count] = sf_open(m_fileOut, SFM_WRITE, &(file_info));
            if(outputfd[count] == NULL)
            {
                char tmp_string[1024]; sprintf(tmp_string, "ERROR: Individual output file %s open error\n", m_fileOut);
                m_log.print_this(tmp_string, 10);
                exit(1);
            }
        }
    }

    ////////////open the delays output files//////////

    //we open the delays output file
    char delays_out[512];
    sprintf(delays_out,"%s%s.del",(*m_config).RESULTPATH.c_str(), (*m_config).SHOWNAME.c_str());
    delfd = fopen(delays_out, "w");
    if(delfd == NULL)
    {
        char tmp_string[1024]; sprintf(tmp_string, "ERROR: output delays file %s open error\n", delays_out);
        m_log.print_this(tmp_string, 10);
        exit(1);
    }

    //we open the second delays output file
    sprintf(delays_out,"%s%s.del2",(*m_config).RESULTPATH.c_str(), (*m_config).SHOWNAME.c_str());
    delfd2 = fopen(delays_out, "w");
    if(delfd2 == NULL)
    {
        char tmp_string[1024]; sprintf(tmp_string, "output delays file %s open error\n", delays_out);
        m_log.print_this(tmp_string, 10);
        exit(1);
    }


    /////////////open the overlap output files////////////

    //we open an overlap detector file where I store the found overlap regions
    if((*m_config).DO_OUTPUT_OVERLAP)
    {
        char overlap_out[512];
        sprintf(overlap_out,"%s%s.ovl",(*m_config).RESULTPATH.c_str(), (*m_config).SHOWNAME.c_str());
        ovlfd = fopen(overlap_out, "w");
        if(ovlfd == NULL)
        {
            char tmp_string[1024]; sprintf(tmp_string, "ERROR: overlaps file %s open error\n", overlap_out);
            m_log.print_this(tmp_string, 10);
            exit(1);
        }
    }

    ////////////open the general info output file//////////
    char info_out[512];
    sprintf(info_out,"%s%s.info",(*m_config).RESULTPATH.c_str(), (*m_config).SHOWNAME.c_str());
    infofd = fopen(info_out, "w");
    if(infofd == NULL)
    {
        char tmp_string[1024]; sprintf(tmp_string, "ERROR: info file %s open error\n", info_out);
        m_log.print_this(tmp_string, 10);
        exit(1);
    }

    ////////////open the weight features output file//////////
    char featw_out[512];
    sprintf(featw_out,"%s%s.weat",(*m_config).RESULTPATH.c_str(), (*m_config).SHOWNAME.c_str());
    featwfd = fopen(featw_out, "w");
    if(featwfd == NULL)
    {
        char tmp_string[1024]; sprintf(tmp_string, "ERROR: info file %s open error\n", featw_out);
        m_log.print_this(tmp_string, 10);
        exit(1);
    }

}

/*!
  Close all used channels and files
  \return Nothing
*/

void FileInOut::Close_Channels()
{
    for(int count=0; count < m_numCh; count++)
    {
        sf_close(inputfd[count]);
        if((*m_config).INDIV_CHANNELS)
        {
            sf_close(outputfd[count]);
        }
    }
    sf_close(outfd[0]);
    if((*m_config).output_second_best_sph)
        sf_close(outfd[1]);
    fclose(delfd);
    fclose(delfd2);
    if((*m_config).DO_OUTPUT_OVERLAP)
        fclose(ovlfd);
    fclose(infofd);
    fclose(featwfd);
}



/*!
  Print the first and second delays (and their xcorr) to the file

  \param finalDelays delays to write to file
  \param finalXcorrValues xcorr values to write to file
  \param UEMGap UEM gap to apply to the file time start
  \return Nothing
*/

void FileInOut::delays_to_file(vector<vector<vector<int> > > & finalDelays, vector<vector<vector<float> > > & finalXcorrValues, float UEMGap)
{
  long frame;
  unsigned int count;
  for(count=0; count<finalDelays[0].size(); count++)
  {
    //print the real frame number
    frame = (*m_config).rate*count+(long)(UEMGap);
    fprintf(delfd,"%ld ->", frame);
    fprintf(delfd2,"%ld ->", frame);
    //for each channel print their delays
    for(int channel_count=0; channel_count<m_numCh; channel_count++)
    {
      fprintf(delfd," %d %f ",finalDelays[channel_count][count][0], finalXcorrValues[channel_count][count][0]);
      fprintf(delfd2," %d %f ",finalDelays[channel_count][count][1], finalXcorrValues[channel_count][count][1]);

    }
    fprintf(delfd,"\n");
    fprintf(delfd2,"\n");
  }
  //we fill up the rest of the delays to match the length of the file by replicating the final ones
  //NOTE that, in addition, we can add some extra delays to ensure the number of delays always exceeds the length of the audio file
  //This is useful for some applications combining delays and acoustic features.
  if((*m_config).DO_DELAYS_PADDING)
  {
      unsigned int full_amount_delays = (unsigned int)((m_frames - UEMGap)/((*m_config).rate*m_sampleRateInMs)) + (*m_config).extra_delays_padding;
      for(count = finalDelays[0].size(); count<full_amount_delays; count++)
      {
        frame = (*m_config).rate*count+(long)(UEMGap);
        fprintf(delfd,"%ld ->", frame);
        fprintf(delfd2,"%ld ->", frame);
        //for each channel print their delays
        for(int channel_count=0; channel_count<m_numCh; channel_count++)
        {
          fprintf(delfd," %d %f ",finalDelays[channel_count][finalDelays[0].size()-1][0], finalXcorrValues[channel_count][finalDelays[0].size()-1][0]);
          fprintf(delfd2," %d %f ",finalDelays[channel_count][finalDelays[0].size()-1][1], finalXcorrValues[channel_count][finalDelays[0].size()-1][1]);

        }
        fprintf(delfd,"\n");
        fprintf(delfd2,"\n");
      }
  }
}


/*!
  Read the acoustic data from the channel files
  \param channel channel number to read from
  \param chanData Where to put the read data
  \param StartSample starting sample to read
  \param numSamples Number of samples to read
*/
int FileInOut::readChanData(int channel, float* chanData, long startSample, int numSamples)
{
    long seeked_frame = sf_seek(inputfd[channel], startSample, SEEK_SET);

    if( seeked_frame == -1)
    {
      char tmp_string[1024]; sprintf(tmp_string,"ERROR: Error trying to read beyond file's reach in channel %d\n", channel);
      m_log.print_this(tmp_string, 10);
      exit(1);
    }
    else if(seeked_frame != startSample)
    {
      char tmp_string[1024]; sprintf(tmp_string,"Error positioning the read pointer for channel %d. Wanted %ld, got %ld\n", channel, startSample, seeked_frame);
      m_log.print_this(tmp_string, 10);
      exit(1);
    }

	return(sf_read_float(inputfd[channel], chanData, numSamples));
}
