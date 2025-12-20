#pragma once

#include <optional>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include "src/usher_graph.hpp"
#include "read.hpp"
#include "util.hpp"

class dataset {
    boost::program_options::variables_map options;
    tbb::global_control init;
public:
    dataset(boost::program_options::variables_map options) :
        options{options}, 
        init{
            tbb::global_control::max_allowed_parallelism,
            options["threads"].as<uint32_t>()
        }
    { }

    uint32_t num_threads() const {
        return options["threads"].as<uint32_t>();
    }
    
    uint32_t max_reads() const {
        return options["max-reads"].as<uint32_t>();
    }
    
    uint32_t min_phred() const {
        return options["min-phred"].as<uint32_t>();
    }
    
    uint32_t clade_idx() const {
        return options["clade-idx"].as<uint32_t>();
    }

    double min_af() const {
        return std::stof(options["min-af"].as<std::string>());
    }
    
    uint32_t min_depth() const {
        return options["min-depth"].as<uint32_t>();
    }
    
    double min_prop() const {
        return std::stof(options["min-prop"].as<std::string>());
    }
    
    std::string wepp_directory() const {
        return options["working-directory"].as<std::string>();
    }
    
    std::string dataset_name() const {
        return options["dataset"].as<std::string>();
    }

    std::string data_directory() const {
        return "./data/" + dataset_name() + "/";
    }

    std::string intermediate_directory() const {
        return "./intermediate/" + dataset_name() + "/";
    }

    std::string results_directory() const {
        return "./results/" + dataset_name() + "/";
    }

    std::string file_prefix() const {
        return options["file-prefix"].as<std::string>();
    }

    std::string ref_path() const {
        return this->data_directory() + options["ref-fasta"].as<std::string>();
    }

    std::string residual_mutations_path() const {
        return this->intermediate_directory() + "residual_mutations.txt";
    }
    
    std::string haplotype_tsv_path() const {
        return this->results_directory() + this->file_prefix() + "_haplotypes.tsv";
    }
    
    std::string haplotype_bam_path() const {
        return this->results_directory() + this->file_prefix() + "_haplotypes.bam";
    }

    const std::vector<int>& masked_sites() const {
        static std::vector<int> cached_mask; 
        if (cached_mask.empty()) {
            std::ifstream inputFile(this->data_directory() + "/mask.bed");
            if (!inputFile.is_open()) {
                // assume no masks
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
    
    std::string checkpoint_path() const {
        return this->intermediate_directory() + this->file_prefix() + "_checkpoint.txt";
    }

    std::string haplotype_read_path() const {
        return this->results_directory() + this->file_prefix() + "_haplotype_reads.csv";
    }

    std::string mutation_reads_path() const {
        return this->results_directory() + this->file_prefix() + "_mutation_reads.csv";
    }
    
    std::string mutation_haplotypes_path() const {
        return this->results_directory() + this->file_prefix() + "_mutation_haplotypes.csv";
    }

    std::string haplotype_proportion_path() const {
        return this->results_directory() + this->file_prefix() + "_haplotype_abundance.csv";
    }

    std::string haplotype_uncertainty_path() const {
        return this->results_directory() + this->file_prefix() + "_haplotype_uncertainty.csv";
    }

    std::string lineage_proportion_path() const {
        return this->results_directory() + this->file_prefix() + "_lineage_abundance.csv";
    }

    std::string haplotype_growth_path() const {
        return this->results_directory() + this->file_prefix() + "_" + this->file_prefix() + "_haplotype_growth.csv";
    }

    const std::string& reference_name() const {
        static std::string ref_header;
        if (ref_header.empty()) {
            std::ifstream fasta_f(this->ref_path());
            if (!fasta_f.is_open())
            {
                std::cerr << "Error: Unable to open file " << ref_path() << std::endl;
                exit(1);
            }
            std::string line;
            std::getline(fasta_f, line);
        
            // Remove leading '>' and only extracting first token
            if (!line.empty() && line[0] == '>') {
                line = line.substr(1);
                // Extract only the first token (accession)
                std::istringstream iss(line);
                iss >> ref_header;
            }
            else {
                std::cerr << "Error: Fasta format NOT correct " << ref_path() << std::endl;
                exit(1);
            }
        }

        return ref_header;
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
                // convert nucleotides to uppercase
		std::transform(temp.begin(), temp.end(), temp.begin(),
			[](unsigned char c){ return std::toupper(c); });
		ref_seq += temp;
            }
            saved.emplace(ref_seq);
        }
        
        return saved.value();
    }

    std::string sam_path() const {
        return this->intermediate_directory() + this->file_prefix() + "_alignment.sam";
    }

    std::string pb_path() const {
        return this->intermediate_directory() + this->file_prefix() + "_reads.pb";
    }
    
    std::vector<raw_read> reads() const;
    std::unordered_map<std::string, std::vector<std::string>> read_reverse_merge() const;

    MAT::Tree mat() const {
        MAT::Tree T;
        T.root = NULL;
        std::string mat_filename = this->data_directory() + this->options["input-mat"].as<std::string>();
        T = MAT::load_mutation_annotated_tree(mat_filename);
        T.uncondense_leaves();
        return T;
    }

    std::vector<std::pair<std::string, std::vector<MAT::Mutation>>>  true_haplotypes() const {
        return read_sample_vcf(this->data_directory() + this->file_prefix() + "_samples.vcf");
    }
};
