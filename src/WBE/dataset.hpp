#pragma once

#include <optional>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include "src/usher_graph.hpp"
#include "read.hpp"
#include "util.hpp"

class dataset {
    boost::program_options::variables_map options;

public:
    dataset(boost::program_options::variables_map options) :
        options{options}
    { }

    uint32_t num_threads() const {
        return options["threads"].as<uint32_t>();
    }

    std::string directory() const {
        return options["output-directory"].as<std::string>() + '/';
    }

    std::string file_prefix() const {
        return options["output-files-prefix"].as<std::string>();
    }

    std::string ref_path() const {
        return this->directory() + options["ref-fasta"].as<std::string>();
    }

    std::string first_checkpoint_path() const {
        return this->directory() + this->file_prefix() + "_first_checkpoint.txt";
    }

    std::string last_checkpoint_path() const {
        return this->directory() + this->file_prefix() + "_last_checkpoint.txt";
    }

    const std::string& reference() const {
        static std::optional<std::string> saved;
        if (!saved) {
            std::ifstream fasta_f(this->ref_path());
            if (!fasta_f.is_open())
            {
                std::cerr << "Error: Unable to open file " << ref_path() << std::endl;
                exit(1);
            }
            std::string ref_header;
            std::getline(fasta_f, ref_header);
            std::string temp;
            std::string ref_seq;
            while (fasta_f)
            {
                std::getline(fasta_f, temp);
                ref_seq += temp;
            }
            saved.emplace(ref_seq);
        }
        
        return saved.value();
    }

    std::string sam_path() const {
        return this->directory() + options["align-sam"].as<std::string>();
    }

    std::string pb_path() const {
        return this->directory() + this->file_prefix() + "_reads.pb";
    }

    std::vector<raw_read> reads() const;

    MAT::Tree mat() const {
        std::string input_mat_filename = this->directory() + this->options["input-mat"].as<std::string>();
        MAT::Tree T = MAT::load_mutation_annotated_tree(input_mat_filename);
        T.uncondense_leaves();
        return T;
    }

    std::vector<std::string> true_haplotypes() const {
        return read_sample_vcf(this->directory() + this->file_prefix() + "_samples.vcf");
    }
};