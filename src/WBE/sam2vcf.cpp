#include "wbe.hpp"

void sam2VCF(po::parsed_options parsed) {
    //main argument for the complex extract command
    po::variables_map vm = parseWBEcommand(parsed);
    std::string dir_prefix = vm["output-directory"].as<std::string>();

    boost::filesystem::path path(dir_prefix);
    if (!boost::filesystem::exists(path)) {
        fprintf(stderr, "Creating output directory.\n\n");
        boost::filesystem::create_directory(dir_prefix);
    }
    path = boost::filesystem::canonical(dir_prefix);
    dir_prefix = path.generic_string();
    dir_prefix += "/";
    std::string vcf_filename_reads = dir_prefix + vm["output-files-prefix"].as<std::string>() + "_reads.vcf";
    std::string freyja_vcf_file = dir_prefix + vm["output-files-prefix"].as<std::string>() + "_reads_freyja.vcf";
    std::string freyja_depth_file = dir_prefix + vm["output-files-prefix"].as<std::string>() + "_reads_freyja.depth";
    std::string ref_fasta = dir_prefix + vm["ref-fasta"].as<std::string>();
    std::string sam_file = dir_prefix + vm["align-sam"].as<std::string>();

    //Reading reference genome
    timer.Start();
    std::ifstream fasta_f(ref_fasta);
    std::string ref_header;
    std::getline(fasta_f, ref_header);
    std::string temp;
    std::string ref_seq;
    while (fasta_f) {
        std::getline(fasta_f, temp);
        ref_seq += temp;
    }

    //Reading sam file
    std::unordered_map<size_t, std::string> read_map;
    std::unordered_map<int, std::vector<std::tuple<std::string, std::string, std::vector<size_t>>>> site_MutReads_map;
    readSAM(sam_file, ref_seq, read_map, site_MutReads_map);

    //Write VCFs
    std::ofstream outfile_vcf(vcf_filename_reads, std::ios::out | std::ios::binary);
    std::ofstream outfile_freyja_vcf(freyja_vcf_file, std::ios::out | std::ios::binary);
    std::ofstream outfile_freyja_depth(freyja_depth_file, std::ios::out | std::ios::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::output> outbuf_vcf, outbuf_freyja_vcf, outbuf_freyja_depth;
    if (vcf_filename_reads.find(".gz\0") != std::string::npos) {
        outbuf_vcf.push(boost::iostreams::gzip_compressor());
    }
    if (freyja_vcf_file.find(".gz\0") != std::string::npos) {
        outbuf_freyja_vcf.push(boost::iostreams::gzip_compressor());
    }
    if (freyja_depth_file.find(".gz\0") != std::string::npos) {
        outbuf_freyja_depth.push(boost::iostreams::gzip_compressor());
    }
    outbuf_vcf.push(outfile_vcf);
    outbuf_freyja_vcf.push(outfile_freyja_vcf);
    outbuf_freyja_depth.push(outfile_freyja_depth);
    std::ostream vcf(&outbuf_vcf);
    std::ostream freyja_vcf(&outbuf_freyja_vcf);
    std::ostream freyja_depth(&outbuf_freyja_depth);
    
    //Write Header
    vcf << "##fileformat=VCFv4.2\n##reference=stdin:hCoV-19/Wuhan/Hu-1/2019|EPI_ISL_402125|2019-12-31\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    freyja_vcf << "##fileformat=VCFv4.2\n##reference=stdin:hCoV-19/Wuhan/Hu-1/2019|EPI_ISL_402125|2019-12-31\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
    for (size_t i = 0; i < read_map.size(); i++)
        vcf << "\t" << read_map[i];

    //Write Mutations
    for (int i = 1; i <= (int)ref_seq.size(); i++) {
        char ref_nuc = ref_seq[i - 1];
        if (site_MutReads_map.find(i) != site_MutReads_map.end()) {
            size_t rd_count = 0;
            bool mut_found = false, no_indel_mut_found = false;;
            //Check if mutation is present at site and create list of ref_nuc possibilities
            for (auto& mr_tuple: site_MutReads_map[i]) {
                rd_count += std::get<2>(mr_tuple).size();
                if (std::get<0>(mr_tuple) != std::get<1>(mr_tuple)) {
                    mut_found = true;
                    if ((std::get<0>(mr_tuple) == std::string(1, ref_nuc)) && (std::get<1>(mr_tuple).size() == 1))
                        no_indel_mut_found = true;
                }
            }
            //Freyja Depth
            freyja_depth <<  "NC_045512v2\t" + std::to_string(i) + "\t" + std::string(1, ref_nuc) + "\t" + std::to_string(rd_count) + "\n";

            if (mut_found) {
                //Freyja VCF
                for (auto& mr_tuple: site_MutReads_map[i]) {
                    if (std::get<0>(mr_tuple) != std::get<1>(mr_tuple))
                        freyja_vcf << "\nNC_045512v2\t" + std::to_string(i) + "\t" +std::get<0>(mr_tuple) + std::to_string(i) + std::get<1>(mr_tuple) + "\t" + std::get<0>(mr_tuple) + "\t" + std::get<1>(mr_tuple) + "\t.\t.\tAF=" + std::to_string((double)std::get<2>(mr_tuple).size() / (double)rd_count);
                }
            }
            if (no_indel_mut_found) {
                //Only accounting substitutions for WEPP
                //VCF -> ID 
                int alt_count = 0;
                vcf << "\nNC_045512v2\t" + std::to_string(i) + "\t";
                for (auto& mr_tuple: site_MutReads_map[i]) {
                    if ((std::get<0>(mr_tuple) != std::get<1>(mr_tuple)) && (std::get<0>(mr_tuple) == std::string(1, ref_nuc)) && (std::get<1>(mr_tuple).size() == 1)) {
                        if (alt_count)
                            vcf << ",";
                        vcf << std::string(1, ref_nuc) + std::to_string(i) + std::get<1>(mr_tuple);
                        alt_count++;
                    }
                }
                vcf << "\t" + std::string(1, ref_nuc) + "\t";
                //VCF -> ALT
                alt_count = 0;
                for (auto& mr_tuple: site_MutReads_map[i]) {
                    if ((std::get<0>(mr_tuple) != std::get<1>(mr_tuple)) && (std::get<0>(mr_tuple) == std::string(1, ref_nuc)) && (std::get<1>(mr_tuple).size() == 1)) {
                        if (alt_count)
                            vcf << ",";
                        vcf << std::get<1>(mr_tuple);
                        alt_count++;
                    }
                }
                vcf << "\t.\t.\t.\t.";
                //VCF -> read matrix
                std::vector<std::pair<int, size_t>> nuc_read_list;
                alt_count = 1;
                for (auto& mr_tuple: site_MutReads_map[i]) {
                    if ((std::get<0>(mr_tuple) != std::get<1>(mr_tuple)) && (std::get<0>(mr_tuple) == std::string(1, ref_nuc)) && (std::get<1>(mr_tuple).size() == 1)) {
                        for (const auto& rd_idx: std::get<2>(mr_tuple))
                            nuc_read_list.emplace_back(std::make_pair(alt_count, rd_idx));
                        alt_count++;
                    }
                }
                //Sort the reads based on read_map
                std::sort(nuc_read_list.begin(), nuc_read_list.end(), compareIdx);
                auto nr_itr = nuc_read_list.begin();
                for (size_t j = 0; j < read_map.size(); j++) {
                    if (nr_itr != nuc_read_list.end()) {
                        if (j == nr_itr->second) {
                            vcf << "\t" + std::to_string(nr_itr->first);
                            nr_itr++;
                        }
                        else
                            vcf << "\t0";
                    }
                    else 
                        vcf << "\t0";
                }
                nuc_read_list.clear();
            }
        }
        else
            freyja_depth <<  "NC_045512v2\t" + std::to_string(i) + "\t" + std::string(1, ref_nuc) + "\t0\n";
    }
    fprintf(stderr, "\nFiles generated in %ld sec \n\n", (timer.Stop() / 1000));
}
