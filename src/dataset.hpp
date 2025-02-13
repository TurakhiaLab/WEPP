#pragma once

#include <optional>
#include <boost/filesystem.hpp>
#include <boost/iostreams/filter/lzma.hpp>
#include <boost/program_options.hpp>

#include "panman/panmanUtils.hpp"
#include "panman_bridge.hpp"
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

    uint32_t max_reads() const
    {
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

    std::string residual_mutations_path() const {
        return this->directory() + "/residual_mutations.txt";
    }
    
    std::string haplotype_tsv_path() const {
        return this->directory() + this->file_prefix() + "_haplotypes.tsv";
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

    std::string variants_path() const {
        return "./src/Freyja/cwap_variants.tsv";
    }

    std::string depth_path() const {
        return "./src/Freyja/cwap_depth.tsv";
    }

    const panmanUtils::Tree& mat() const {
        static std::optional<panmanUtils::Tree> tree;
        if (!tree.has_value()) {
            std::string input_mat_filename = this->directory() + this->options["input-mat"].as<std::string>();
            std::ifstream fin(input_mat_filename);
            boost::iostreams::filtering_streambuf< boost::iostreams::input> inPMATBuffer;
            inPMATBuffer.push(boost::iostreams::lzma_decompressor());
            inPMATBuffer.push(fin);
            std::istream inputStream(&inPMATBuffer);
            panmanUtils::TreeGroup g(inputStream);
            assert(g.trees.size() == 1);
            tree.emplace(g.trees[0]);
        }
        return tree.value();
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
        // CONSENSUS Sequence is actually reference sequence in our PanMAT
        static std::optional<std::string> saved;
        if (!saved) {
            const panmanUtils::Tree& mat = this->mat();
            coord_converter converter{mat};
            std::string ref = converter.reference;
            saved.emplace(std::move(ref));
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

    std::vector<std::string> true_haplotypes() const {
        return read_sample_vcf(this->directory() + this->file_prefix() + "_samples.vcf");
    }
};