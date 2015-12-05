#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <tuple>
#include <exception>
#include <algorithm>
#include <string.h>
#include "read_param.hpp"

//This function reads the key, value pairs from the config file into a map. 
//The argument is only used to check the keys are what we expected.
std::map<std::string, std::string> SpbConfigParser::get_parameter_map(std::map<std::string, ValueTuple> expected_keys)
{
  std::fstream config;
  config.open(m_config.c_str(), std::fstream::in);
  if ( ! config.is_open() ) {
      throw std::runtime_error("Can't open config in file: "+ m_config);
  }
    
  //Store a name here when it has already been found in the file
  std::map<std::string, std::string> found;

  /* Read the lines */
  while ( config.good() ) {
    std::string line;
    std::getline(config, line);
    std::pair<std::string, std::string> parsed_line = parse_line(line);
    if (!config.good() )
        break;
    if (parsed_line.first.empty())
        continue;
    //Check not a duplicate
    std::map<std::string, std::string>::iterator fit = found.find(parsed_line.first);
    if (fit != found.end()){
        throw std::runtime_error("Key: "+parsed_line.first+" found twice in "+m_config+" with values "+fit->second+" and "+parsed_line.second);
    }
    //Check that we expected this key
    std::map<std::string, ValueTuple>::iterator dit = expected_keys.find(parsed_line.first);
    if (dit == expected_keys.end()){
        throw std::runtime_error(m_config+" contained unexpected key: "+parsed_line.first);
    }
    //Store the name and value in found
    found[parsed_line.first] = parsed_line.second;
  }
  config.close();
  return found;
}

int SpbConfigParser::parameter_parser(std::map<std::string, ValueTuple > datastore)
{
  std::map<std::string, std::string> found = get_parameter_map(datastore);
  //Now we have read the key, value pairs into a map, store them where we want them
  for(std::map<std::string, ValueTuple>::iterator dit = datastore.begin(); dit != datastore.end(); ++dit) {
    //Get the string we want. Use found unless it wasn't, in which case use the default
    std::string value = std::get<2>(dit->second);
    std::map<std::string, std::string>::iterator fit = found.find(dit->first);
    if (fit != found.end() ) {
        value = fit->second;
    }
    else if (m_verbose) {
        std::cout<<"Using default value: "<<value<<" for key: "<<dit->first<<std::endl;
    }
    if (value.empty()) {
        throw std::runtime_error("No value for key: "+dit->first);
    }
    //Switch on the type of the key
    TypeEnum thistype = std::get<1>(dit->second);
    if(thistype == IntType) {
        * (int *) std::get<0>(dit->second) = atoi(value.c_str());
    }
    else if (thistype == StringType) {
        * (std::string *) std::get<0>(dit->second) = value;
    }
    else if (thistype == FloatType) {
        * (double *) std::get<0>(dit->second) = atof(value.c_str());
    }
    else {
        throw std::runtime_error("Unrecognised type for "+dit->first);
    }
  }
  return 0;
}

//Parse a single line into a key-value pair
std::pair<std::string, std::string> SpbConfigParser::parse_line(std::string linein)
{
    std::string key, value;
    //Strip comments
    size_t comment_pos = linein.find_first_of(m_comments);
    linein = linein.substr(0, comment_pos);
    //Strip leading whitespace
    size_t white_pos = linein.find_first_not_of(std::string(" \t"));
    if (white_pos == std::string::npos)
        return std::make_pair(key, value);
    linein = linein.substr(white_pos);
    //Find key/value separator
    size_t separator_pos = linein.find_first_of(m_separators);
    key = linein.substr(0, separator_pos);
    //Set value
    if (separator_pos != std::string::npos){
        size_t value_pos = linein.find_first_not_of(m_separators,separator_pos);
        if (value_pos != std::string::npos){
            value = linein.substr(value_pos);
        }
    }
    //Take out whitespace
    key.erase(std::remove(key.begin(), key.end(), ' '), key.end());
    key.erase(std::remove(key.begin(), key.end(), '\t'), key.end());
    value.erase(std::remove(value.begin(), value.end(), '\t'), value.end());
    value.erase(std::remove(value.begin(), value.end(), ' '), value.end());
    return std::make_pair(key, value);
}
