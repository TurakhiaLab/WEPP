#include <time.h>

#include <array>
#include <boost/program_options.hpp>

#include "util.hpp"
#include "pipeline.hpp"
#include "sam2pb.hpp"
#include "analyze_peaks.hpp"

int main (int argc, char** argv) {
    namespace po = boost::program_options;

    po::options_description global("Command options");
    global.add_options()
    ("command", po::value<std::string>(), "Command to execute. Valid options are annotate, mask, extract, uncertainty, and summary.")
    ("subargs", po::value<std::vector<std::string> >(), "Command-specific arguments.");
    po::positional_options_description pos;
    pos.add("command",1 ).add("subargs", -1);
    std::string cmd;
    po::variables_map vm;
    po::parsed_options parsed = po::command_line_parser(argc, argv).options(global).positional(pos).allow_unregistered().run();
    //this help string shows up over and over, lets just define it once
    std::string cnames[] = {"COMMAND","initial_filter","post_filter", "sam2PB"};
    std::string chelp[] = {
        "DESCRIPTION\n\n",
        "Detects peaks on selected lineages\n\n",
        "Converts given sam to pb for running the pipeline\n\n"
    };
    try {
        po::store(parsed, vm);
        cmd = vm["command"].as<std::string>();
    } catch (...) { //not sure this is the best way to catch it when matUtils is called with no positional arguments.
        fprintf(stderr, "\nNo command selected. Help follows:\n\n");
        for (size_t i = 0; i < std::size(cnames); ++i) {
            fprintf(stderr, "%-15s\t%s", cnames[i].c_str(), chelp[i].c_str());
        }
        //0 when no command is selected because that's what passes tests.
        exit(0);
    }
    if (cmd == "initial_filter") { 
        dataset ds{parseWBEcommand(parsed)};
        detect_peaks(ds, true);
    }
    else if (cmd == "post_filter") { 
        dataset ds{parseWBEcommand(parsed)};
        detect_peaks(ds, false);
    } else if (cmd == "sam2PB") {
        dataset ds{parseWBEcommand(parsed)};
        sam2PB(ds);
    } else if (cmd == "analyzePeaks") {
        dataset ds{parseWBEcommand(parsed)};
        analyze_peaks(ds);
    } else if (cmd == "help") {
        fprintf(stderr, "\n");
        for (size_t i = 0; i < std::size(cnames); ++i) {
            fprintf(stderr, "%-15s\t%s", cnames[i].c_str(), chelp[i].c_str());
        }
        exit(0);
    } else {
        fprintf(stderr, "\nInvalid command. Help follows:\n\n");
        for (size_t i = 0; i < std::size(cnames); ++i) {
            fprintf(stderr, "%-15s\t%s", cnames[i].c_str(), chelp[i].c_str());
        }
        exit(1);
    }

    return 0;
}