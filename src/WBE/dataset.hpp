#pragma once

#include <optional>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include "src/usher_graph.hpp"
#include "read.hpp"
#include "util.hpp"

class dataset {
    boost::program_options::variables_map options;
    tbb::task_scheduler_init init;
public:
    dataset(boost::program_options::variables_map options) :
        options{options}, init{(int) options["threads"].as<uint32_t>()}
    { }

    uint32_t num_threads() const {
        return options["threads"].as<uint32_t>();
    }
    
    uint32_t max_reads() const {
        return options["max-reads"].as<uint32_t>();
    }

    std::string directory() const {
        return options["output-directory"].as<std::string>() + '/';
    }

    std::string comparison_directory() const {
        return options["comparison-directory"].as<std::string>() + '/';
    }

    std::string file_prefix() const {
        return options["output-files-prefix"].as<std::string>();
    }

    std::string comparison_file_prefix() const {
        return options["comparison-files-prefix"].as<std::string>();
    }

    std::string ref_path() const {
        return this->directory() + options["ref-fasta"].as<std::string>();
    }

    std::string residual_mutations_path() const {
        return this->directory() + "/residual_mutations.txt";
    }
    
    std::string haplotype_sam_path() const {
        return this->directory() + this->file_prefix() + "_haplotypes.sam";
    }
    
    std::string haplotype_bam_path() const {
        return this->directory() + this->file_prefix() + "_haplotypes.bam";
    }

    const std::vector<int>& masked_sites() const {
        static std::vector<int> cached_mask; 
        if (cached_mask.empty()) {
            std::ifstream inputFile("./mask.bed");
            if (!inputFile.is_open()) {
                std::cerr << "Error: Unable to open mask.bed file." << std::endl;
                return cached_mask; 
            }

            std::string line;
            while (std::getline(inputFile, line)) {
                std::istringstream lineStream(line);
                std::string col1, col2;
                int col3;

                // Parse the line into columns
                if (lineStream >> col1 >> col2 >> col3) {
                    cached_mask.push_back(col3);
                }
            }
            inputFile.close();
        }
        return cached_mask;
    }
    
    std::string first_checkpoint_path() const {
        return this->directory() + this->file_prefix() + "_first_checkpoint.txt";
    }

    std::string haplotype_read_path() const {
        return this->directory() + this->file_prefix() + "_haplotype_reads.csv";
    }

    std::string mutation_reads_path() const {
        return this->directory() + this->file_prefix() + "_mutation_reads.csv";
    }
    
    std::string mutation_haplotypes_path() const {
        return this->directory() + this->file_prefix() + "_mutation_haplotypes.csv";
    }

    std::string haplotype_proportion_path() const {
        return this->directory() + this->file_prefix() + "_haplotype_abundance.csv";
    }

    std::string comparison_haplotype_proportion_path() const {
        return this->comparison_directory() + this->comparison_file_prefix() + "_haplotype_abundance.csv";
    }

    std::string haplotype_growth_path() const {
        return this->directory() + this->file_prefix() + "_" + this->comparison_file_prefix() + "_haplotype_growth.csv";
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
    
    std::string comparison_pb_path() const {
        return this->comparison_directory() + this->comparison_file_prefix() + "_reads.pb";
    }

    std::vector<raw_read> reads() const;
    std::unordered_map<std::string, std::vector<std::string>> read_reverse_merge() const;

    MAT::Tree mat() const {
        MAT::Tree T;
        T.root = NULL;
        std::string mat_filename = this->directory() + this->options["input-mat"].as<std::string>();
        T = MAT::load_mutation_annotated_tree(mat_filename);
        T.uncondense_leaves();
        return T;
    }

    MAT::Tree cmp_mat() const {
       MAT::Tree T;
       T.root = NULL;

       std::string mat_filename = this->comparison_directory() + this->options["cmp-mat"].as<std::string>();
       T = MAT::load_mutation_annotated_tree(mat_filename);
       T.uncondense_leaves();
       return T;
    }

    std::vector<std::pair<std::string, std::vector<MAT::Mutation>>>  true_haplotypes() const {
        return read_sample_vcf(this->directory() + this->file_prefix() + "_samples.vcf");
    }
};