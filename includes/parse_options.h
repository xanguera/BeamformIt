#include <iostream>
#include <string>
#include <cstring>
#include <map>
#include <stdio.h>
#include <stdlib.h>
#include "global.h"

using namespace std;

//This class holds the parameters we will use in the system
//first they are created by add_param and then filled using command line or config file
class Params
{
public:
    bool set_value(string value);

    //parameter values
    string m_longID; //long ID
    string m_shortID; //short ID
    void* m_value; //where the value goes to
    string m_type; //type of the data the value goes to
    bool m_mandatory; //true when it is mandatory
    bool m_externallySet; //true when this parameter has been externally set
    string m_default; //default value
    string m_param_desc; //description of this parameter
};

class parse_options
{
public:
  //constructor setting up the config file
  parse_options(Configuration* config)
  {
    m_config = config;
  }

  int add_param(string longID, string shortID, string type, void* destination, bool mandatory, string defaultValue, string param_desc);

  void print_help();

  //parses the input params, including those from an input file
  void parse_params(int argc, char *argv[]);

  void processConfigFile(string configFile);
  bool processCommandLine(int argc, char *argv[]);
  
private:

  Configuration* m_config;
  map<string, Params*> m_params; //stores the input parameters and points to the values (key is always the long version)
  map<string, std::string > m_options; //mapping for input options to the option name used to handle short and long versions. All go to the long version
};
