#include "haplotype_pruning.hpp"

po::variables_map parse_haplotype_pruning_command(po::parsed_options parsed) {
    po::variables_map vm;
    po::options_description conv_desc("haplotype_pruning options");
    conv_desc.add_options()
    ("input-mat,i", po::value<std::string>()->required(),
     "Input mutation-annotated tree file [REQUIRED]")
    ("output-mat,w", po::value<std::string>()->default_value(""),
     "Write the selected tree as a new protobuf to the target file.")
    ("output-directory,o", po::value<std::string>()->default_value("./"),
     "Write output files to the target directory. Default is current directory.")
    ("samples,s", po::value<std::string>()->default_value(""),
     "Select samples by explicitly naming them. One per line")
    ("mut-distance,d", po::value<int>()->default_value(-1),
     "Select samples which have mutation distance less than or equal to d from samples.")
    ("help,h", "Print help messages");
    // Collect all the unrecognized options from the first pass. This will include the
    // (positional) command name, so we need to erase that.
    std::vector<std::string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
    opts.erase(opts.begin());

    // Run the parser, with try/catch for help
    try {
        po::store(po::command_line_parser(opts)
                  .options(conv_desc)
                  .run(), vm);
        po::notify(vm);
    } catch(std::exception &e) {
        std::cerr << conv_desc << std::endl;
        // Return with error code 1 unless the user specifies help
        if (vm.count("help"))
            exit(0);
        else
            exit(1);
    }
    return vm;
}


void haplotype_pruning(po::parsed_options parsed) {
    po::variables_map vm = parse_haplotype_pruning_command(parsed);
    std::string dir_prefix = vm["output-directory"].as<std::string>();
    std::string input_mat_filename = dir_prefix + vm["input-mat"].as<std::string>();
    std::string output_mat_filename = dir_prefix + vm["output-mat"].as<std::string>();
    std::string input_samples_file = dir_prefix + vm["samples"].as<std::string>();
    int mut_dist = vm["mut-distance"].as<int>();

    // Load input MAT and uncondense tree
    MAT::Tree T;
    T = MAT::load_mutation_annotated_tree(input_mat_filename);
    T.uncondense_leaves();

    //Getting sample selection from file
    std::vector<std::string> samples;
    if (input_samples_file != "") {
        samples = read_sample_names(input_samples_file);
    }

    //Get remaining samples
    samples = samples_outside_mut_dist(T, samples, mut_dist);

    //Filter the input based on the samples
    MAT::Tree subtree = filter_master(T, samples, false, true);
    auto n = T.get_node("node_102181");
    int d = mutation_distance(T, n, T.get_node("England/QEUH-9CA801/2020|OA983298.1|2020-09-10"));
    printf("Orig distance Node_102181: %d\n",d);


    subtree.condense_leaves();
    MAT::save_mutation_annotated_tree(subtree, output_mat_filename);
}


//Function to remove haplotypes within certain distance from given haplotypes
std::vector<std::string> samples_outside_mut_dist(const MAT::Tree &T, std::vector<std::string> samples_to_check, int mut_dist_thresh) {
    std::vector<std::string> good_samples;
    for (auto n: T.depth_first_expansion()) {
        if (n->is_leaf()) {
            bool remove = false;
            for (auto sample: samples_to_check) {
                //Remove leaf node within mut_dist threshold of sample
                int m_dist = mutation_distance(T, T.get_node(sample), n);
                if (m_dist <= mut_dist_thresh) {
                    remove = true;
                    break;
                }
            }
            if (!remove)
                good_samples.emplace_back(n->identifier);
        }
    }
    return good_samples;
}