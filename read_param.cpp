#include <fstream>
#include <string>
#include <map>
#include <algorithm>
#include "read_param.hpp"

std::vector<std::string> SpbConfigParser::GetRemainingKeys()
{
    std::vector<std::string> keys;
    //C++11!
    for(auto pair: optdict)
        keys.push_back(pair.first);
    return keys;
}

//This function reads the key, value pairs from the config file into a map. 
//The argument is only used to check the keys are what we expected.
std::map<std::string, std::string> SpbConfigParser::get_parameter_map()
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
    //Store the name and value in found
    found[parsed_line.first] = parsed_line.second;
  }
  config.close();
  return found;
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

template<> int SpbConfigParser::convert_value<int>(std::string value)
{
    return atoi(value.c_str());
}

template<> double SpbConfigParser::convert_value<double>(std::string value)
{
    return atof(value.c_str());
}

//For strings just return them as is
template<> std::string SpbConfigParser::convert_value<std::string>(std::string value)
{
    return value;
}

template<> bool SpbConfigParser::convert_value<bool>(std::string value)
{
    return (bool) atoi(value.c_str());
}
