#pragma once

#include <optional>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include "src/usher_graph.hpp"
#include "read.hpp"

class dataset {
    std::string directory;
    std::string file_prefix;
    std::optional<std::string> fasta_path;
    std::string sam_file;

public:
    dataset(po::parsed_options options) :
        options{options}
    { }

    std::string ref_name() const {

    }

    std::string sam_name() const {

    }

    std::string pb_name() const {

    }

    std::vector<read> reads() const {
        
    }

    MAT::Tree mat() const {

    }

    std::vector<std::string> true_haplotypes() {

    }
};