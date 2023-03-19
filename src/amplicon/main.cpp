#include "amplicon.hpp"
#include "src/mutation_annotated_tree.hpp"
#include "tbb/concurrent_unordered_set.h"
#include <algorithm>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <cassert>
#include <random>
#include <stdio.h>
#include <string>
#include <time.h>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace po = boost::program_options;

static void print_mutation_in_range(const std::string &sample_name,
                                    const MAT::Node *sample_node,
                                    int start_range, int end_range) {
    std::unordered_map<int, MAT::Mutation> mutations;
    while (sample_node) {
        for (const auto &mut : sample_node->mutations) {
            if (mut.position >= start_range && mut.position <= end_range) {
                mutations.emplace(mut.position, mut);
            }
        }
        sample_node = sample_node->parent;
    }
    for (const auto &mut : mutations) {
        if (mut.second.ref_nuc != mut.second.mut_nuc) {
            printf("DEBUG:\t%s\t%d\t%c\n", sample_name.c_str(), mut.first,
                   MAT::get_nuc(mut.second.mut_nuc));
        }
    }
}

int main(int argc, char **argv) {
    Timer timer;

    std::string outdir;

    po::options_description desc("optimize options");
    desc.add_options()("input-mat,i", po::value<std::string>()->required(),
                       "Input mutation-annotated tree file [REQUIRED].")(
        "outdir,d", po::value<std::string>(&outdir)->default_value("."),
        "Output directory to dump output and log files [DEFAULT uses current "
        "directory]")(
        "vcf,v", po::value<std::string>()->required(),
        "Input VCF file (in uncompressed or gzip-compressed .gz format) "
        "[REQUIRED]")("fasta,f", po::value<std::string>()->required(),
                      "Input aligned fasta alignment file) "
                      "[REQUIRED]");

    po::options_description all_options;
    all_options.add(desc);
    po::positional_options_description p;
    po::variables_map vm;
    try {
        po::store(po::command_line_parser(argc, argv)
                      .options(all_options)
                      .positional(p)
                      .run(),
                  vm);
        po::notify(vm);
    } catch (std::exception &e) {
        std::cerr << desc << std::endl;
        // Return with error code 1 unless
        // the user specifies help
        if (vm.count("help"))
            exit(0);
        else
            exit(1);
    }
    std::string input_mat_filename = vm["input-mat"].as<std::string>();
    std::string input_vcf_filename = vm["vcf"].as<std::string>();
    std::string input_fasta_filename = vm["fasta"].as<std::string>();

    timer.Start();
    fprintf(stdout, "Loading input MAT file = %s\n",
            input_mat_filename.c_str());
    fprintf(stdout, "Loading input FASTA file = %s\n",
            input_fasta_filename.c_str());

    // Load input MAT and uncondense tree
    MAT::Tree T;
    T = MAT::load_mutation_annotated_tree(input_mat_filename);
    T.uncondense_leaves();
    fprintf(stderr, "Completed in %ld msec \n", timer.Stop());

    // All input amplicon sequences to place, not found already on
    // input tree
    std::vector<Missing_Sample> missing_samples;

    // Read VCF of missing amplicons and extract new amplicon sequences to place
    MAT::read_vcf(&T, input_vcf_filename, missing_samples,
                  false); // create_new_mat = false currently

    auto num_samples = missing_samples.size();
    // Check that there are actually missing samples from VCF to add to tree
    if (num_samples <= 0) {
        printf("No new samples found in input VCF to add to input MAT. "
               "Exiting program.\n");
        return 0;
    }
    printf("Found %lu missing amplicon samples.\n", num_samples);

    // Amplicon samples to place (debugging, can comment out)
    //    printf("Printing missing amplicon samples place on input MAT\n");
    //    for (auto &sample : missing_samples) {
    //        printf("%s\n", sample.name.c_str());
    //    }

    // Creating output directories
    boost::filesystem::path path(outdir);
    if (!boost::filesystem::exists(path)) {
        boost::filesystem::create_directory(path);
    }
    path = boost::filesystem::canonical(outdir);
    outdir = path.generic_string();
    std::vector<FILE*> annotation_files;
    for (auto annot_idx=0; annot_idx<T.get_num_annotations(); annot_idx++) {
        annotation_files.push_back(fopen((outdir+"/clade"+std::to_string(annot_idx)+".txt").c_str(), "w"));
    }

    static tbb::affinity_partitioner ap;

    // <start, end> coordinate for each amplicon sequence in input vcf
    // (currently assumes not sorting amplicon samples)
    std::vector<std::tuple<int, int>> samples_start_end_coordinates;

    // Get start and ending coordinates of each amplicon aligned to reference
    // from aligned multifasta file (input given as -f arg)
    get_coordinates(input_fasta_filename, samples_start_end_coordinates);

    /* DEBUGGING for amplicon start/end coordinate correctness from alignment
    for (auto &tup : sample_start_end_coordinates) {
      printf("START: %d\t", std::get<0>(tup));
      printf("END: %d\n", std::get<1>(tup));
    }
    */
    printf("All amplicon alignment coordinates retrieved from %s\n",
           input_fasta_filename.c_str());

    //  Holds best parsimony scores and number of parsimony-optimal placement
    //  for each new amplicon placed on tree
    std::vector<int> best_parsimony_scores;
    std::vector<size_t> num_best_placements;
    auto num_annotations = T.get_num_annotations();

    printf("Starting placement of each amplicon sample iteratively on tree.\n");
    // Iterate over all the missing amplicon samples
    for (size_t s = 0; s < missing_samples.size(); ++s) {

        // For now, assuming indexes are 0,1,2... this is the order that
        // the amplicon samples will be placed sequentially
        // on the tree, in the order they are retrieved from vcf (according to
        // read_vcf) and placed into missing_samples vector. Update this later
        // with another vector of indexes into the missing_samples vector and
        // the samples_start_end_coordinates vector if we want to sort amplicon
        // samples before placement.

        //  Sort the missing sample mutations by position
        std::sort(missing_samples[s].mutations.begin(),
                  missing_samples[s].mutations.end());

        auto bfs = T.breadth_first_expansion();
        size_t total_nodes = bfs.size();

        std::vector<std::vector<MAT::Mutation>> node_excess_mutations(
            total_nodes);
        std::vector<std::vector<MAT::Mutation>> imputed_mutations(total_nodes);

        // Parsimony score of placing each sample at each node of the tree (DFS
        // order)
        std::vector<int> node_set_difference(total_nodes);

        size_t best_node_num_leaves = 0;
        int best_set_difference = 1e9;
        std::vector<bool> node_has_unique(total_nodes);
        size_t best_j = 0;
        bool best_node_has_unique = false;
        size_t best_distance = 1e9;
        std::vector<size_t> best_j_vec;
        size_t num_best = 1;
        MAT::Node *best_node = T.root;
        best_j_vec.emplace_back(0);
        std::vector<size_t> node_distance(total_nodes, 0);

        // Searching for most parsimonious placements.
        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, total_nodes),
            [&](tbb::blocked_range<size_t> r) {
                for (size_t k = r.begin(); k < r.end(); ++k) {

                    node_has_unique[k] = false;
                    mapper2_input inp;
                    inp.T = &T;
                    inp.node = bfs[k];
                    inp.missing_sample_mutations =
                        &missing_samples[s].mutations;
                    inp.excess_mutations = &node_excess_mutations[k];
                    inp.imputed_mutations = &imputed_mutations[k];
                    inp.best_node_num_leaves = &best_node_num_leaves;
                    inp.best_set_difference = &best_set_difference;
                    inp.best_node = &best_node;
                    inp.best_j = &best_j;
                    inp.num_best = &num_best;
                    inp.j = k;
                    inp.has_unique = &best_node_has_unique;

                    inp.set_difference = &node_set_difference[k];

                    inp.distance = node_distance[k];
                    inp.best_distance = &best_distance;

                    inp.best_j_vec = &best_j_vec;
                    inp.node_has_unique = &(node_has_unique);

                    mapper2_body(inp, true);
                }
            },
            ap);

        fprintf(stdout, "Currently placing amplicon: %s\n",
                missing_samples[s].name.c_str());

        // Get starting and ending coordinates for masking amplicon sample
        int start = std::get<0>(samples_start_end_coordinates[s]);
        int end = std::get<1>(samples_start_end_coordinates[s]);

        int min_parsimony = best_set_difference;
        auto sample = missing_samples[s].name;
        num_best = 0;
        std::vector<size_t> best_placements;
        for (size_t k = 0; k < total_nodes; k++) {
            size_t num_mut = 0;

            // Is placement as sibling
            if (bfs[k]->is_leaf() || node_has_unique[k]) {
                std::vector<MAT::Mutation> common_mut, l1_mut, l2_mut;
                std::vector<MAT::Mutation> curr_l1_mut;

                for (auto m1 : bfs[k]->mutations) {
                    MAT::Mutation m = m1.copy();
                    curr_l1_mut.emplace_back(m);
                }

                // Compute l2_mut
                for (auto m1 : node_excess_mutations[k]) {
                    if ((m1.position < start) || (m1.position > end)) {
                        MAT::Mutation m = m1.copy();
                        common_mut.emplace_back(m);
                        continue;
                    }
                    bool found = false;
                    for (auto m2 : curr_l1_mut) {
                        if (m1.is_masked()) {
                            break;
                        }
                        if ((m1.position == m2.position) &&
                            (m1.mut_nuc == m2.mut_nuc)) {
                            found = true;
                            MAT::Mutation m = m1.copy();
                            common_mut.emplace_back(m);
                            break;
                        }
                    }
                    if (!found) {
                        MAT::Mutation m = m1.copy();
                        l2_mut.emplace_back(m);
                    }
                }
                // Compute l1_mut
                for (auto m1 : curr_l1_mut) {
                    bool found = false;
                    for (auto m2 : common_mut) {
                        if (m1.is_masked()) {
                            break;
                        }
                        if (m1.position == m2.position) {
                            if (m1.mut_nuc == m2.mut_nuc) {
                                found = true;
                                break;
                            }
                        }
                    }
                    if (!found) {
                        MAT::Mutation m = m1.copy();
                        l1_mut.emplace_back(m);
                    }
                }

                num_mut = l2_mut.size();

                if (num_mut < min_parsimony) {
                    best_j = k;
                    best_placements.clear();
                    best_placements.push_back(k);
                    min_parsimony = num_mut;
                    num_best = 1;
                } else if (num_mut == min_parsimony) {
                    best_placements.push_back(k);
                    num_best++;
                }
            }
            // Else placement as child
            else {
                std::vector<MAT::Mutation> node_mut;

                std::vector<MAT::Mutation> curr_l1_mut;

                for (auto m1 : bfs[k]->mutations) {
                    MAT::Mutation m = m1.copy();
                    curr_l1_mut.emplace_back(m);
                }

                for (auto m1 : node_excess_mutations[k]) {
                    bool found = false;
                    if ((m1.position < start) || (m1.position > end)) {
                        continue;
                    }
                    for (auto m2 : curr_l1_mut) {
                        if (m1.is_masked()) {
                            break;
                        }
                        if ((m1.position == m2.position) &&
                            (m1.mut_nuc == m2.mut_nuc)) {
                            found = true;
                            break;
                        }
                    }
                    if (!found) {
                        MAT::Mutation m = m1.copy();
                        node_mut.emplace_back(m);
                    }
                }
                num_mut = node_mut.size();

                if (num_mut < min_parsimony) {
                    best_j = k;
                    min_parsimony = num_mut;
                    best_placements.clear();
                    best_placements.push_back(k);
                    num_best = 1;
                } else if (num_mut == min_parsimony) {
                    best_placements.push_back(k);
                    num_best++;
                }
            }
        }
        assert(num_best == best_placements.size());
        printf("Best placement parsimony score: %i\n", min_parsimony);
        printf("Number of equally parsimonious placements: %zu\n", num_best);

        best_node = bfs[best_j];

        //std::vector<std::vector<std::string>> clade_assignments;
        std::vector<std::unordered_map<std::string, int>> clade_assignments(num_annotations);
        for (size_t c = 0; c < num_annotations; c++) {
            auto& this_annotation_class=clade_assignments[c];
            for (size_t k = 0; k < best_placements.size(); k++) {
                bool include_self =
                    !bfs[best_placements[k]]->is_leaf() && !node_has_unique[k];
                auto clade_assignment =
                    T.get_clade_assignment(bfs[best_placements[k]], c, include_self);
                auto res=this_annotation_class.emplace(clade_assignment,1);
                if (!res.second) {
                    res.first->second++;
                }
            }
        }
        for (size_t k = 0; k < num_annotations; k++) {
            auto annotations_file=annotation_files[k];
            fprintf(annotations_file, "%s\t", sample.c_str());
            std::vector<std::pair<std::string, int>> clade_counts(clade_assignments[k].begin(),clade_assignments[k].end());
            std::sort(clade_counts.begin(),clade_counts.end(),[](const std::pair<std::string, int>& first, const std::pair<std::string, int>& second){
                if(first.second>second.second){
                    return true;
                }else if (first.second==second.second) {
                    return first.first>second.first;
                }
                return false;
            });
            for(int clade_idx=0;clade_idx<clade_counts.size();clade_idx++){
                const auto& clade_to_print=clade_counts[clade_idx];
                fprintf(annotations_file, "%s(%i/%zu,%0.2f%%)",
                                clade_to_print.first.c_str(), clade_to_print.second,
                                num_best,100*(float)clade_to_print.second/(float)num_best);
                if (clade_idx!=clade_counts.size()-1) {
                    fputc(',', annotations_file);
                }
            }
            
            fprintf(annotations_file, "\n");
            fflush(annotations_file);
        }
        // Is placement as sibling
        if (best_node->is_leaf() || node_has_unique[best_j]) {
            std::string nid = T.new_internal_node_id();
            T.create_node(nid, best_node->parent->identifier);
            T.create_node(sample, nid);
            T.move_node(best_node->identifier, nid);
            // common_mut stores mutations common to the
            // best node branch and the sample, l1_mut
            // stores mutations unique to best node branch
            // and l2_mut stores mutations unique to the
            // sample not in best node branch
            std::vector<MAT::Mutation> common_mut, l1_mut, l2_mut;
            std::vector<MAT::Mutation> curr_l1_mut;

            // Compute current best node branch mutations
            for (auto m1 : best_node->mutations) {
                MAT::Mutation m = m1.copy();
                curr_l1_mut.emplace_back(m);
            }
            // Clear mutations on the best node branch which
            // will be later replaced by l1_mut
            best_node->clear_mutations();

            // Compute l2_mut
            for (auto m1 : node_excess_mutations[best_j]) {
                if ((m1.position < start) || (m1.position > end)) {
                    MAT::Mutation m = m1.copy();
                    common_mut.emplace_back(m);
                    continue;
                }
                bool found = false;
                for (auto m2 : curr_l1_mut) {
                    if (m1.is_masked()) {
                        break;
                    }
                    if ((m1.position == m2.position) &&
                        (m1.mut_nuc == m2.mut_nuc)) {
                        found = true;
                        MAT::Mutation m = m1.copy();
                        common_mut.emplace_back(m);
                        break;
                    }
                }
                if (!found) {
                    MAT::Mutation m = m1.copy();
                    l2_mut.emplace_back(m);
                }
            }
            // Compute l1_mut
            for (auto m1 : curr_l1_mut) {
                bool found = false;
                for (auto m2 : common_mut) {
                    if (m1.is_masked()) {
                        break;
                    }
                    if (m1.position == m2.position) {
                        if (m1.mut_nuc == m2.mut_nuc) {
                            found = true;
                            break;
                        }
                    }
                }
                if (!found) {
                    MAT::Mutation m = m1.copy();
                    l1_mut.emplace_back(m);
                }
            }

            // Add mutations to new node using common_mut
            for (auto m : common_mut) {
                T.get_node(nid)->add_mutation(m);
            }
            // Add mutations to best node using l1_mut
            for (auto m : l1_mut) {
                T.get_node(best_node->identifier)->add_mutation(m);
            }
            // Add new sample mutations using l2_mut
            for (auto m : l2_mut) {
                T.get_node(sample)->add_mutation(m);
            }
        }
        // Else placement as child
        else {
            T.create_node(sample, best_node->identifier);
            MAT::Node *node = T.get_node(sample);
            std::vector<MAT::Mutation> node_mut;

            std::vector<MAT::Mutation> curr_l1_mut;

            for (auto m1 : best_node->mutations) {
                MAT::Mutation m = m1.copy();
                curr_l1_mut.emplace_back(m);
            }

            for (auto m1 : node_excess_mutations[best_j]) {
                bool found = false;
                if ((m1.position < start) || (m1.position > end)) {
                    continue;
                }
                for (auto m2 : curr_l1_mut) {
                    if (m1.is_masked()) {
                        break;
                    }
                    if ((m1.position == m2.position) &&
                        (m1.mut_nuc == m2.mut_nuc)) {
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    MAT::Mutation m = m1.copy();
                    node_mut.emplace_back(m);
                }
            }
            for (auto m : node_mut) {
                node->add_mutation(m);
            }
        }
#ifndef NDEBUG
        print_mutation_in_range(sample, T.get_node(sample), start, end);
#endif
    }

    bool retain_original_branch_len = false;
    auto final_tree_filename = outdir + "/final-tree.nh";
    fprintf(stderr, "Writing final tree to file %s \n",
            final_tree_filename.c_str());
    auto parsimony_score = T.get_parsimony_score();
    /* Print parsimony score for tree
       fprintf(stderr, "The parsimony score for this tree is: %zu \n",
       parsimony_score);
    */

    std::ofstream final_tree_file(final_tree_filename.c_str(),
                                  std::ofstream::out);
    std::stringstream newick_ss;
    write_newick_string(newick_ss, T, T.root, true, true,
                        retain_original_branch_len);
    final_tree_file << newick_ss.rdbuf();
    final_tree_file.close();

    fprintf(stderr, "Completed in %ld msec \n\n", timer.Stop());
}

