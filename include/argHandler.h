//
//  argHandler.h
//  EvaporationSim
//
//  Created by Andreas von Haaren on 21/05/2020.
//  Copyright © 2020 Andreas von Haaren. All rights reserved.
//

#ifndef argHandler_h
#define argHandler_h

#include <filesystem>                   // filesystem
namespace fs = std::filesystem;

#include <map>                          // dictionary functionality
#include <iostream>                     // input + output
#include "dsmc.h"                       // DSMC implementation


namespace args {
    /// Enum class of validation codes for the configuration file
    enum class config_validity {
        valid,                      ///< enum value for a valid config
        fileNotExist,               ///< enum value for a non-existent config file
        incorrectFormat,            ///< enum value for an incorrect file format
        invalidOtherReason          ///< enum value for any other reason
    };
    /// Function to check the validity of a configuration file
    /// @param configFile a const reference to a path to the configuration file
    /// @returns a config_validity value
    inline config_validity checkConfigValidity(const fs::path& configFile) {
        if (!fs::exists(configFile))
            return config_validity::fileNotExist;
        if (configFile.extension() != ".json")
            return config_validity::incorrectFormat;
        if (!fs::is_regular_file(configFile))
            return config_validity::invalidOtherReason;
        return config_validity::valid;
    }

    /// Enum class of validation codes for the output folder
    enum class output_validity {
        valid,                      ///< enum value for a valid folder
        dirNotExist,                ///< enum value for non-existent folder
        notADirectory,              ///< enum value if the supplied path is not a folder
        dirNotEmpty                 ///< enum value for a non-empty folder
    };
    /// Function to check the validity of an output folder
    /// @param out_path a const reference to a path to the output folder
    /// @returns an output_validity value
    inline output_validity checkOutputDirValidity(const fs::path& out_path) {
        if (not fs::exists(out_path))
            return output_validity::dirNotExist;
        if (not fs::is_directory(out_path))
            return output_validity::notADirectory;
        if (not fs::is_empty(out_path))
            return output_validity::dirNotEmpty;
        return output_validity::valid;
    }
    
    /// Handler for the arguments to main()
    /// @discussion Checks the configuration file and output folder and applies any options that can be passed.
    /// @param argc the number of arguments + 1
    /// @param argv a const char* array with the arguments to main()
    /// @param options a map from string to bool reference that contains the possible options that can be passed to main()
    /// @warning Exits with code 1 if any of the necessary arguments are invalid
    inline void handle(int argc, const char* argv[], std::map<std::string, bool&> options = {}) {
        // Lambda to print the Usage message
        auto printUsage = [&argv]() {
            static const std::string usage_msg = "\nUsage: \n\n" +
            std::string(argv[0]) +
            " {options} [path_to_config.json] [path_to_output]\n"
            "Options:\n"
            "    -ow: overwrite existing files without asking for confirmation\n"
            "    -animate=[fps]: output necessary files for animating the cloud motion\n"
            "                    [fps] = the number of frames per second (integer), default is 60\n"
            "    -seed: seed random generators to get randomised results\n"
            "\n";
            std::cout << usage_msg;
        };
        
        std::string overwriteOption = "-ow";
        std::string animateOption = "-animate";
        
        // Fill a vector with the arguments (not options)
        std::vector<std::string> arguments;
        for (int i = 1; i < argc; ++i) {
            // options are marked with a '-'
            if (argv[i][0] == '-') {
                std::string opt = argv[i];
                if (opt.find(animateOption) == std::string::npos) {
                    try {
                        options.at(opt) = true;
                    } catch (std::out_of_range& oor) {
                        std::cout << "'" << opt.substr(1, opt.size() - 1) << "' is not a valid option.\n";
                    }
                } else {
                    if (opt.substr(0, animateOption.size()) != animateOption) {
                        std::cout << "'" << opt.substr(1, opt.size() - 1) << "' is not a valid option.\n";
                    } else {
                        options.at(animateOption) = true;
                        if (opt.size() > 8) {
                            try {
                                DSMC::frametime = 1./std::stod(opt.substr(9, opt.size() - 1));
                            } catch (const std::exception& e) {
                                std::cout << "The animate option has to be in the format '" << animateOption << "[=fps]' where fps is an integer value.\n";
                                exit(1);
                            }
                        } else {
                            DSMC::frametime = 1./60;
                        }
                    }
                }
            } else {
                arguments.push_back(argv[i]);
            }
        }
        
        // there are 2 necessary arguments, the configuration file and the output folder
        if (arguments.size() == 2) {
            // Configuration file is expected first, output folder second
            std::string first_arg = arguments[0], second_arg = arguments[1];
            // the configuration file is only checked superficially. The contents are only checked when it is read
            auto check = checkConfigValidity(first_arg);
            using cval = config_validity;
            switch (check) {
                case cval::valid:
                    break;
                case cval::fileNotExist:
                    std::cout << "\n" << first_arg << " does not exist\n";
                    printUsage();
                    exit(1);
                case cval::incorrectFormat:
                    std::cout << "\nThe specified configuration file is not a json file.\n";
                    printUsage();
                    exit(1);
                case cval::invalidOtherReason:
                    std::cout << "\nSomething went wrong.\n\n";
                    exit(1);
                default:
                    break; // This should never be reached.
            }
            
            // check the output folder
            auto dirCheck = checkOutputDirValidity(second_arg);
            using oval = output_validity;
            std::string confirm;
            switch (dirCheck) {
                case oval::valid:
                    break;
                case oval::dirNotExist:
                    std::cout << "\n" << second_arg << " does not exist.\n\n";
                    exit(404);
                case oval::notADirectory:
                    std::cout << "\n" << second_arg << " is not a directory.";
                    printUsage();
                    exit(1);
                case oval::dirNotEmpty:
                    // if the folder has any contents and the overwrite option is not set,
                    // the user will be asked if they want to overwrite the contents
                    if (options.at(overwriteOption))
                        break;
                    else {
                        std::cout << "\nContents of " << second_arg << " will be overwritten. Proceed [Y/N]?   ";
                        std::getline(std::cin, confirm);
                        if (confirm == "Y" || confirm == "y")
                            break;
                        else {
                            std::cout << "Exiting.\n\n";
                            exit(1);
                        }
                    }
                default:
                    break; // This should never be reached.
            }
            
            // If everything is okay, set the respective static members of DSMC
            DSMC::configFile = fs::path(first_arg);
            DSMC::outputFolder = fs::path(second_arg);
        }
        else {
            printUsage();
            exit(1);
        }
    }

    
}


#endif /* argHandler_h */
