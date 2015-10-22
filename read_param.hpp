#ifndef READ_PARAM_H
#define READ_PARAM_H

#include <tuple>
#include <string>
#include <map>

enum TypeEnum
{
    StringType, IntType, FloatType
};

typedef std::tuple<void *, enum TypeEnum, std::string>  ValueTuple;

/* 
 * I will become the nth person to reimplement a simple ini-file parser in C++ 
 * The existing wheels are not round enough!
 * This will set all keys not found in the config file to their default value.
 * If a key is found twice, it is an error.
 * Finding unexpected keys is also an error
 * Input values:
 *      separators: string containing all characters to use as name: value separator characters. 
 *      Old Gadget format has just a space, we also want =
 *      comments: string containing comment characters. Old Gadget uses %, we also want #
 *      std::map<std::string name, std::tuple<void * ptr, enum ValueType type, std::string default> > - 
 *      associates config file names with locations, types, and default values.
 *      Search for std::string name in the file, and set the location pointed to by ptr
 *      */
class SpbConfigParser
{
    public:
        SpbConfigParser(const char * config_filename, const std::string& separators=" =\t", const std::string& comments="%#", bool verbose=false): m_verbose(verbose), m_config(config_filename), m_separators(separators), m_comments(comments)
        {
        }

        //This function reads the key, value pairs from the config file into a map. 
        //The argument is only used to check the keys are what we expected.
        std::map<std::string, std::string> get_parameter_map(std::map<std::string, ValueTuple> expected_keys);
        //Parses the whole file and stores all the values in the addresses pointed to by the pointers in datastore 
        int parameter_parser(std::map<std::string, ValueTuple > datastore);
        //Parse a single line into a key-value pair
        std::pair<std::string, std::string> parse_line(std::string linein);
    private:
        const bool m_verbose;
        const std::string m_config;
        const std::string m_separators;
        const std::string m_comments;
};
#endif
