#ifndef READ_PARAM_H
#define READ_PARAM_H

#include <string>
#include <map>
#include <exception>
#include <iostream>
#include <vector>

/* 
 * I will become the nth person to reimplement a simple ini-file parser in C++ 
 * The existing wheels are not round enough!
 * If a key is found twice, it is an error.
 * Input values:
 *      separators: string containing all characters to use as name: value separator characters. 
 *      Old Gadget format has just a space, we also want =
 *      comments: string containing comment characters. Old Gadget uses %, we also want #
 *      */
class SpbConfigParser
{
    public:
        SpbConfigParser(const std::string& config_filename, const std::string& separators=" =\t", const std::string& comments="%#", bool verbose=false): m_verbose(verbose), m_config(config_filename), m_separators(separators), m_comments(comments)
        {
            optdict = get_parameter_map();
        }
        //Get a value from the parser dictionary, removing it afterwards.
        //In this version, raise an exception if key is not found
        template<typename T> T PopValue(std::string key)
        {
            std::string value;
            std::map<std::string, std::string>::iterator fit = optdict.find(key);
            if (fit != optdict.end() ) {
                value = fit->second;
            }
            else {
                throw std::domain_error("Key: "+key+" not found");
            }
            if (value.empty()) {
                throw std::invalid_argument("No value for key: "+key);
            }
            return convert_value<T>(value);
        }

        //In this version, use the default if key is not found
        template<typename T> T PopValue(std::string key, T deflt)
        {
            try {
                return PopValue<T>(key);
            }
            catch(std::domain_error){
                if (m_verbose)
                    std::cout<<"Using default value: "<<deflt<<" for key: "<<key<<std::endl;
                return deflt;
            }
        }

        //Get all keys that were found but have not yet been called with PopValue
        std::vector<std::string> GetRemainingKeys();
    private:
        //This function reads the key, value pairs from the config file into a map. 
        std::map<std::string, std::string> get_parameter_map();
        //Parse a single line into a key-value pair
        std::pair<std::string, std::string> parse_line(std::string linein);
        //Convert a string into a more useful type
        template<typename T> T convert_value(std::string value)
        {
            //This is the default: types we recognise get specialisations.
            throw std::runtime_error("Unrecognised type");
        }
        const bool m_verbose;
        const std::string m_config;
        const std::string m_separators;
        const std::string m_comments;
        std::map<std::string, std::string> optdict;
};

//Specialisations defined in read_param.cpp.
template<> int SpbConfigParser::convert_value<int>(std::string value);
template<> double SpbConfigParser::convert_value<double>(std::string value);
template<> std::string SpbConfigParser::convert_value<std::string>(std::string value);
template<> bool SpbConfigParser::convert_value<bool>(std::string value);

#endif
