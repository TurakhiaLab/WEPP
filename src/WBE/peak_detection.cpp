#include "wbe.hpp"

constexpr bool USE_FAST = true;
constexpr int TOP_N = 25;
constexpr int TOP_N_EXPANDED_FACTOR = 4;
constexpr int MAX_NEIGHBORS = 100;
constexpr int MAX_PEAK_PEAK_MUTATION = 1;
constexpr int MAX_PEAK_NONPEAK_MUTATION = 4;
constexpr int MAX_CACHED_MULTIPLICITY_SIZE = 1024;
constexpr int PARALLEL_READS_THRESHOLD = 256;

void detectPeaks(po::parsed_options parsed) {
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
    std::string input_mat_filename = dir_prefix + vm["input-mat"].as<std::string>();
    std::string prior_lineages = vm["prior-lineages"].as<std::string>();
    std::string vcf_filename_samples = dir_prefix + vm["output-files-prefix"].as<std::string>() + "_samples.vcf";
    std::string proto_reads = dir_prefix + vm["output-files-prefix"].as<std::string>() + "_reads.pb";
    std::string hap_csv_filename = dir_prefix + vm["output-files-prefix"].as<std::string>() + "_haplotype_abundance.csv";
    std::string freyja_lineage_csv_filename = dir_prefix + vm["output-files-prefix"].as<std::string>() + "_freyja_results.csv";
    std::string barcode_file = dir_prefix + vm["output-files-prefix"].as<std::string>() + "_barcode.csv";
    std::string condensed_nodes_csv = dir_prefix + vm["output-files-prefix"].as<std::string>() + "_condensed_nodes.csv";
    std::string read_mutation_depth_vcf = dir_prefix + vm["output-files-prefix"].as<std::string>() + "_read_data.vcf";
    std::string ref_fasta = dir_prefix + vm["ref-fasta"].as<std::string>();
    uint32_t num_threads = vm["threads"].as<uint32_t>();

    fprintf(stderr, "\nNum Cores: %d\n\n", num_threads);
    tbb::task_scheduler_init init(num_threads);
    
    timer.Start();
    //Get prior lineages
    std::vector<std::string> selected_lineages; 
    std::stringstream lin_str(prior_lineages);
    std::string str;
    while (std::getline(lin_str,str,',')) {
        selected_lineages.emplace_back(str);
    }
    //Loading reference genome
    std::ifstream fasta_f(ref_fasta);
    if (!fasta_f.is_open()) {
        std::cerr << "Error: Unable to open file " << ref_fasta << std::endl;
        exit(1); 
    }
    std::string ref_header;
    std::getline(fasta_f, ref_header);
    std::string temp;
    std::string ref_seq;
    while (fasta_f) {
        std::getline(fasta_f, temp);
        ref_seq += temp;
    }
    fprintf(stderr, "\nLoading input MAT files %s \n", input_mat_filename.c_str());

    // Load input MAT and uncondense tree
    MAT::Tree T = MAT::load_mutation_annotated_tree(input_mat_filename);
    T.uncondense_leaves();
    fprintf(stderr, "Completed in %ld sec \n\n", (timer.Stop() / 1000));

    //Read samples.vcf to check how close are peaks to samples 
    std::vector<std::string> vcf_samples;
    readSampleVCF(vcf_samples, vcf_filename_samples);
    //Get the input reads data
    std::unordered_map<size_t, struct read_info*> read_map;
    std::unordered_map<std::string, std::vector<std::string>> reverse_merge;
    load_reads_from_proto(proto_reads, read_map, reverse_merge);
    
    //Get haplotype abundances and condensed node names
    timer.Start();
    std::unordered_map<std::string, double> hap_abun_map, freyja_lineage_abun_map;
    std::unordered_map<std::string, std::vector<std::string>> condensed_nodeNames_map;
    std::vector<MAT::Node*> curr_peak_nodes;
    readCSV(hap_abun_map, hap_csv_filename);
    readCSV(condensed_nodeNames_map, condensed_nodes_csv);
    
    //Get Freyja Lineages
    readCSV(freyja_lineage_abun_map, freyja_lineage_csv_filename);

    //CREATE condensed tree for neighbor lineage search
    int neighbor_dist_thresh = 4;
    MAT::Tree T_condensed;
    std::unordered_map<MAT::Node*, std::vector<MAT::Node*>> condensed_node_mappings;
    createCondensedTree(T.root, read_map, condensed_node_mappings, T_condensed);
    
    //EXTRACT lineages from haplotype names
    std::string check_string = "CONDENSED";
    for (const auto& hap_abun: hap_abun_map) {
        if (hap_abun.first.find(check_string) != std::string::npos) {
            auto node_names_list = condensed_nodeNames_map[hap_abun.first];
            for (int i = 0; i < (int)node_names_list.size(); i++) {
                auto node_name = node_names_list[i];
                size_t last_underscore = node_name.find_last_of('_');
                std::string real_node_name = node_name.substr(0, last_underscore);
                while (T.get_node(real_node_name) == NULL) {
                    last_underscore = real_node_name.find_last_of('_');
                    real_node_name = real_node_name.substr(0, last_underscore);
                }
                if (!i)
                    curr_peak_nodes.emplace_back(T_condensed.get_node(real_node_name));
                std::string curr_lineage = node_name.substr(last_underscore + 1);
                if (std::find(selected_lineages.begin(), selected_lineages.end(), curr_lineage) == selected_lineages.end())
                    selected_lineages.emplace_back(curr_lineage);
            }

        }
        else {
            size_t last_underscore = hap_abun.first.find_last_of('_');
            std::string real_node_name = hap_abun.first.substr(0, last_underscore);
            while (T.get_node(real_node_name) == NULL) {
                last_underscore = real_node_name.find_last_of('_');
                real_node_name = real_node_name.substr(0, last_underscore);
            }
            curr_peak_nodes.emplace_back(T_condensed.get_node(real_node_name));
            std::string curr_lineage = hap_abun.first.substr(last_underscore + 1);
            if (std::find(selected_lineages.begin(), selected_lineages.end(), curr_lineage) == selected_lineages.end())
                selected_lineages.emplace_back(curr_lineage);
        }
    }

    //ADD Neighboring LINEAGES
    addNeighborLineages(T_condensed, T, condensed_node_mappings, curr_peak_nodes, selected_lineages, neighbor_dist_thresh);
    curr_peak_nodes.clear();
    condensed_node_mappings.clear();
    MAT::clear_tree(T_condensed);

    //ADD Freyja lineages
    for (const auto& lin_abun: freyja_lineage_abun_map) {
        auto& curr_lineage = lin_abun.first;
        if (std::find(selected_lineages.begin(), selected_lineages.end(), curr_lineage) == selected_lineages.end()) {
            selected_lineages.emplace_back(curr_lineage);
            printf("Freya lineage: %s\n", curr_lineage.c_str());
        }
    }
    fprintf(stderr, "Added Lineages in %ld sec\n\n", (timer.Stop() / 1000));
    
    //CREATE new tree containg only selected_lineages
    tbb::concurrent_hash_map<MAT::Node*, double> node_score_map;
    MAT::Tree T_lin;
    createLineageTree(T.root, selected_lineages, T_lin);
    selected_lineages.clear();

    analyzeReads(T, T_lin, ref_seq, read_map, node_score_map, vcf_samples, barcode_file, read_mutation_depth_vcf, condensed_nodes_csv);
    
    //RUN peaks_filtering
    std::string command = "python src/WBE/peaks_filtering.py " + vm["output-files-prefix"].as<std::string>() + " " + std::string(dir_prefix) + " 100";
    int result = std::system(command.c_str());
    if (result)
        fprintf(stderr, "\nCannot run peak_filtering.py\n");
}

struct AuxNode {
    double score;
    int leaf_count;
    bool mapped;
    /* is there possibly a single nonmapped (or neighor mapped)?
       node in the subtree? */
    bool subtree_alive;
    bool is_leaf;

    /* children and mutations  */
    AuxNode* parent;
    /* muts from root to here */
    std::vector<MAT::Mutation> stack_muts;
    std::string id;
    MAT::Node* condensed_source;
    std::vector<AuxNode*> children;

    std::vector<int> mutations(std::vector<MAT::Mutation> const& comp, int min_pos, int max_pos) {
        int i = 0, j = 0;
        std::vector<int> muts;

        int last_i = stack_muts.size();
        while (i < (int) stack_muts.size() && stack_muts[i].position < min_pos) ++i;
        while (last_i > 0 && stack_muts[last_i - 1].position > max_pos) --last_i;

        while (i < last_i || j < (int) comp.size()) {
            if (i == last_i) {                
                muts.push_back(comp[j].position);
                ++j;
            }
            else if (stack_muts[i].position < min_pos) {
                ++i;
            }
            else if (stack_muts[i].position > max_pos) {
                return muts;
            }
            else if (j == (int) comp.size()) {
                muts.push_back(stack_muts[i].position);
                ++i;
            }
            else if (stack_muts[i].position < comp[j].position) {
                muts.push_back(stack_muts[i].position);
                ++i;
            }
            else if (stack_muts[i].position > comp[j].position) {
                muts.push_back(comp[j].position);
                ++j;
            }
            else if (stack_muts[i].position == comp[j].position && stack_muts[i].mut_nuc != comp[j].mut_nuc) {
                muts.push_back(comp[j].position);
                ++i; ++j;
            }
            else {
                ++i; ++j;
            }
        }

        return muts;
    }

    int mutation_distance(std::vector<MAT::Mutation> const& comp, int min_pos, int max_pos) {
        return mutations(comp, min_pos, max_pos).size();
    }

    /* N not counted */
    int mutation_distance(struct read_info* other) {
        return mutation_distance(other->mutations, other->start, other->end);
    }

    /* N is counted */
    int mutation_distance(AuxNode * other) {
        return mutation_distance(other->stack_muts, 0, INT_MAX);
    }
};

struct ScoreComparator {
    bool operator() (AuxNode* const& left, AuxNode* const& right) const {
        if (abs(left->score - right->score) > 1e-6) {
            return left->score > right->score;
        }
        else if (left->leaf_count != right->leaf_count) {
            return left->leaf_count > right->leaf_count;
        }
        else {
            return left->id > right->id;
        }
    }
};

struct MutationComparator {
    bool operator() (AuxNode* const& left, AuxNode* const& right) const {
        if (left->stack_muts.size() != right->stack_muts.size()) {
            return left->stack_muts.size() < right->stack_muts.size();
        }  
        for (size_t i = 0; i < left->stack_muts.size(); ++i) {
            if (left->stack_muts[i].position != right->stack_muts[i].position) {
                return left->stack_muts[i].position < right->stack_muts[i].position;
            }
            else if (left->stack_muts[i].mut_nuc != right->stack_muts[i].mut_nuc) {
                return left->stack_muts[i].mut_nuc < right->stack_muts[i].mut_nuc;
            }
        }

        return false;
    }
};

/* slightly larger than O(NM) */
class AuxManager
{
    std::vector<AuxNode> arena;
    const std::vector<read_info *>& reads;
    std::vector<int> max_parismony;
    std::vector<int> parsimony_multiplicity;
    
    /* if epp[i].size() != multiplicity[i], that means not cached */
    /* since there's too many. Also, values are arena indices */
    std::vector<std::set<int>> epp_positions_cache;
    
    /* mappable nodes based on their scores */
    std::set<AuxNode*, ScoreComparator> current_nodes;

    /* indices of unmapped reads */
    std::set<int> remaining_reads;

    /* map of read_index to built */
    std::map<int, AuxNode*> correspondence;

    std::vector<MAT::Node*> peaks, neighbors;
    std::set<AuxNode*, MutationComparator> selected_peaks, selected_neighbors;

    /* can only map if one or mutation positions was removed */
    /* (even if total size increases) */
    int mutation_reductions(std::vector<int> const& parent, std::vector<int> const& child) {
        int reductions = 0;
        for (size_t i = 0, j = 0; i < parent.size(); ++i) {
            while (j < child.size() && child[j] < parent[i]) {
                ++j;
            }
            if (j == child.size() || child[j] != parent[i]) {
                ++reductions;
            }
            else {
                ++j;
            }
        }

        return reductions;
    }

    /* max_val corresponds to max parismony */
    /* maps single read to entire tree, caching aliveness */
    bool single_read_tree(AuxNode* curr, std::vector<int> const & parent_locations, struct read_info *read, std::set<int> &max_indices, int &max_val) {
        std::vector<int> const my_locations = curr->mutations(read->mutations, read->start, read->end);

        /* no need to worry about dead subtrees */
        if (!curr->subtree_alive) {
            return false;
        }
        else if (!curr->mapped) {
            /* map self */
            // this basically ensures semantics are the exact same
            // as the previous algorithm
            int parsimony, reductions = mutation_reductions(parent_locations, my_locations);
            if (reductions) {
                parsimony = parent_locations.size() - reductions;
            }
            else {
                parsimony = my_locations.size();
            }

            if (!curr->is_leaf || reductions) {
                if (parsimony < max_val)
                {
                    max_val = parsimony;
                    max_indices = {(int)(curr - &arena[0])};
                }
                else if (parsimony == max_val)
                {
                    max_indices.insert(curr - &arena[0]);
                }
            }
        }

        bool any_alive = !curr->mapped;
        for (AuxNode *child: curr->children) {
            any_alive |= single_read_tree(child, my_locations, read, max_indices, max_val);
        }

        return curr->subtree_alive = any_alive;
    }

    /* map all reads to all nodes, maintaining caches if not too large */
    void cartesian_map() {
        using my_mutex_t = tbb::queuing_mutex;
        my_mutex_t my_mutex;

        std::vector<int> empty_mutation_list;
        
        epp_positions_cache.resize(reads.size());
        max_parismony.resize(reads.size());
        parsimony_multiplicity.resize(reads.size());

        tbb::parallel_for(tbb::blocked_range<size_t>(0, reads.size()), 
            [&](tbb::blocked_range<size_t> k) {
            for (size_t r = k.begin(); r != k.end(); ++r) {
                std::set<int> max_indices;
                int max_val = INT32_MAX; /* really want minimum value */

                single_read_tree(&arena[0], empty_mutation_list, reads[r], max_indices, max_val);

                double const delta = (double) reads[r]->degree / std::log2(1 + max_indices.size());
                {
                    my_mutex_t::scoped_lock my_lock{my_mutex};
                    for (int arena_index: max_indices) {
                        arena[arena_index].score += delta;
                    }
                }

                this->max_parismony[r] = max_val;
                this->parsimony_multiplicity[r] = max_indices.size();
                if (max_indices.size() <= MAX_CACHED_MULTIPLICITY_SIZE) {
                    this->epp_positions_cache[r] = std::move(max_indices);
                }
            }
        });

        for (size_t i = 0; i < arena.size(); ++i) {
            this->current_nodes.insert(&arena[i]);
        }
    }

    /* find all current reads that map to this node */
    // a bit annoying that there's duplicate code
    std::vector<int> find_correspondents(AuxNode *node) {
        int const arena_index = node - &arena[0];
        std::vector<int> correspondents;

        using my_mutex_t = tbb::queuing_mutex;
        my_mutex_t my_mutex;

        // if large enough, parallelize
        if (remaining_reads.size() >= PARALLEL_READS_THRESHOLD) {
            tbb::parallel_for(tbb::blocked_range<size_t>(0, reads.size()), [&](tbb::blocked_range<size_t> r) {
                for (size_t read = r.begin(); read != r.end(); ++read) {
                    if (this->remaining_reads.find(read) != this->remaining_reads.end()) {
                        if (epp_positions_cache[read].find(arena_index) != epp_positions_cache[read].end()) {
                            my_mutex_t::scoped_lock lock{my_mutex};
                            correspondents.push_back(read); 
                        }
                        else if (epp_positions_cache[read].size() == (size_t) parsimony_multiplicity[read]) {
                            /* cached and not found */
                            continue;
                        }
                        else {
                            struct read_info *r = this->reads[read];
                            std::vector<int> parent = node->parent ? node->parent->mutations(r->mutations, r->start, r->end) : std::vector<int>();
                            std::vector<int> ours = node->mutations(r->mutations, r->start, r->end);

                            int parsimony, reductions = mutation_reductions(parent, ours);
                            if (reductions) {
                                parsimony = parent.size() - reductions;
                            }
                            else {
                                parsimony = ours.size();
                            }
                            
                            if ((!node->is_leaf || reductions) && parsimony == max_parismony[read]) {
                                my_mutex_t::scoped_lock lock{my_mutex};
                                correspondents.push_back(read);
                            }
                        }
                    }
                }
            });
            std::sort(correspondents.begin(), correspondents.end());
        }
        else {
            for (int read: this->remaining_reads) {
                if (epp_positions_cache[read].find(arena_index) != epp_positions_cache[read].end()) {
                    correspondents.push_back(read); 
                }
                else if (epp_positions_cache[read].size() == (size_t) parsimony_multiplicity[read]) {
                    /* cached and not found */
                    continue;
                }
                else {
                    struct read_info *r = this->reads[read];
                    std::vector<int> parent = node->parent ? node->parent->mutations(r->mutations, r->start, r->end) : std::vector<int>();
                    std::vector<int> ours = node->mutations(r->mutations, r->start, r->end);

                    int parsimony, reductions = mutation_reductions(parent, ours);
                    if (reductions) {
                        parsimony = parent.size() - reductions;
                    }
                    else {
                        parsimony = ours.size();
                    }
                    
                    if ((!node->is_leaf || reductions) && parsimony == max_parismony[read]) {
                        correspondents.push_back(read);
                    }
                }
            }
        }
        return correspondents;
    }

    /* removes a singular read's contribution from current tree */
    void remove_reads_effects(int index) {
        std::set<int> recalcuated;
        std::set<int> *epps;
        if (this->epp_positions_cache[index].size() == (size_t) parsimony_multiplicity[index]) {
            epps = &this->epp_positions_cache[index];
        }
        else {
            /* recalculate on (almost) entire tree, generally the most optimal */
            int max_val = INT32_MAX;
            std::vector<int> empty;
            single_read_tree(&arena[0], empty, reads[index], recalcuated, max_val);
            epps = &recalcuated;
        }

        /* use ORIGINAL size */
        double const delta = (double) reads[index]->degree / std::log2(1 + parsimony_multiplicity[index]);
        for (int arena_index: (*epps)) {
            AuxNode *node = &arena[arena_index];
            if (node->mapped) {
                continue;
            }
            this->current_nodes.erase(node);
            node->score -= delta;
            if (node->score > 1e-6) {
                this->current_nodes.insert(node);
            }
        }
        this->remaining_reads.erase(index);
    }

    void singular_step(AuxNode* node) {
        assert(!node->mapped);
        
        /* 1. map remaining reads onto this node (finding correspondents) */
        std::vector<int> correspondents = find_correspondents(node);

        /* 2. map all said reads onto entire remaining set, adjusting deltas */
        /* 3. mark correspondence */ 
        for (int correspondent: correspondents) {
            this->remove_reads_effects(correspondent);
            correspondence[correspondent] = node;
        }
    }

    /* can both aux nodes be sent in the same batch? */
    bool valid_two_tops(AuxNode *a, AuxNode *b) {
         return a->mutation_distance(b) > MAX_PEAK_PEAK_MUTATION;
    }

    bool dfs_possible_neighbors(AuxNode* pivot, AuxNode *curr, std::set<AuxNode*, ScoreComparator>& s, int radius) {
        if (!curr || !curr->subtree_alive) {
            return false;
        }
        else if (pivot->mutation_distance(curr) > radius) {
            return curr->subtree_alive;
        }

        bool alive = false;
        if (!curr->mapped) {
            s.insert(curr);
            alive = true;
        }

        for (AuxNode* child: curr->children) {
            alive |= dfs_possible_neighbors(pivot, child, s, radius);
        }

        return curr->subtree_alive = alive;
    }

    void possible_neighbors(AuxNode* pivot, std::set<AuxNode*, ScoreComparator>& s, int radius) {
        AuxNode *curr = pivot;
        while (curr->parent && pivot->mutation_distance(curr->parent) <= radius) {
            curr = curr->parent;
        }

        dfs_possible_neighbors(pivot, curr, s, radius);
    }

    /* clears neighbors (and self) from map */
    void clear_neighbors(std::vector<AuxNode*> &nodes) {
        for (AuxNode* const& n: nodes) {
            current_nodes.erase(n);
            selected_peaks.insert(n);
            peaks.emplace_back(n->condensed_source);
        }

        std::vector<AuxNode*> added;
        /* candidates */
        for (AuxNode *pivot: nodes) {
            std::set<AuxNode*, ScoreComparator> multisource_radius;
            possible_neighbors(pivot, multisource_radius, MAX_PEAK_NONPEAK_MUTATION);

            int i = 0;
            for (AuxNode * node : multisource_radius)
            {
                bool const found = selected_peaks.find(node) != selected_peaks.end() || selected_neighbors.find(node) != selected_neighbors.end();
                if (!found) {
                    node->mapped = true;
                    current_nodes.erase(node);
                    added.emplace_back(node);
                    neighbors.emplace_back(node->condensed_source);
                    if (++i == MAX_NEIGHBORS)
                    {
                        break;
                    }
                }
            }
        }

        for (AuxNode *pivot: nodes) {
            std::set<AuxNode*, ScoreComparator> multisource_radius;
            possible_neighbors(pivot, multisource_radius, MAX_PEAK_PEAK_MUTATION);

            for (AuxNode * node : multisource_radius)
            {
                node->mapped = true;
                current_nodes.erase(node);
            }
        }

        selected_neighbors.insert(added.begin(), added.end());
    }

    /* matches currently top_n nodes with their reads */
    /* returns if finished or not */
    bool step(int top_n) {
        assert(!current_nodes.empty() && !remaining_reads.empty());

        /* to get exactly same behavior as before, have a single considering set */
        std::vector<AuxNode*> consideration;

        /* get top_n (under some special constraints) */ 
        auto it = current_nodes.begin();
        double const min_score = (*it)->score;

        if (min_score < 1e-6) {
            return true;
        }

        for (int i = 0; 
                i < top_n * TOP_N_EXPANDED_FACTOR && 
                it != current_nodes.end() && 
                abs((*it)->score - min_score) < 1e-6 &&
                consideration.size() < (size_t) top_n; 
            ++i, ++it) 
        {
            bool valid = true;
            for (AuxNode *old: consideration) {
                if (!valid_two_tops(old, *it)) {
                    valid = false;
                }
            }

            if (valid) {
                consideration.push_back(*it);
                (*it)->mapped = true;
                printf("%.9f id: %s\n", (*it)->score, (*it)->id.c_str());
            }
        }

        this->clear_neighbors(consideration);
        
        for (AuxNode *&node: consideration) {
            this->singular_step(node);
        }

        return remaining_reads.empty() || current_nodes.empty();
    }

    AuxNode* from_mat_tree(AuxNode* parent, MAT::Node* node, std::unordered_map<MAT::Node*, std::vector<MAT::Node*>>& condensed_node_mappings) {
        this->arena.emplace_back();
        AuxNode *ret = &this->arena.back();
        ret->parent = parent;
        ret->mapped = false;
        ret->subtree_alive = true;
        ret->is_leaf = condensed_node_mappings.at(node).front()->is_leaf(); 
        ret->score = 0;
        ret->stack_muts = {};
        ret->id = node->identifier;
        ret->condensed_source = node;
        ret->leaf_count = getNumLeaves(condensed_node_mappings, node);
        if (parent) {
            for (const MAT::Mutation &mut : parent->stack_muts)
            {
                /* only need to compare to our muts since parent's are unique */
                bool valid = true;
                for (const MAT::Mutation &comp: node->mutations) {
                    if (comp.position == mut.position) {
                        valid = false;
                        break;
                    }
                }
                
                if (valid) ret->stack_muts.push_back(mut);
            }
        }
        /* maybe rewinded mutation back to original */
        for (size_t i = 0; i < node->mutations.size(); ++i) {
            if (node->mutations[i].ref_nuc != node->mutations[i].mut_nuc) {
                ret->stack_muts.push_back(node->mutations[i]);
            }
        }
        std::sort(ret->stack_muts.begin(), ret->stack_muts.end());

        for (MAT::Node* child : node->children) {
            AuxNode *curr = from_mat_tree(ret, child, condensed_node_mappings); 
            ret->children.emplace_back(curr);
        }

        return ret;
    }

    int mat_tree_size(MAT::Node *node) {
        if (!node) return 0;
        int ret = 1;
        for (MAT::Node* child: node->children) {
            ret += mat_tree_size(child);
        }
        return ret;
    }

public:
    /* preprocessing */
    AuxManager(const MAT::Tree &src, std::unordered_map<MAT::Node*, std::vector<MAT::Node*>>& condensed_node_mappings, std::vector<read_info*> const& reads) : reads{reads} {
        /* create auxiliary tree, ensure no reallocations */
        this->arena.reserve(mat_tree_size(src.root));
        this->from_mat_tree(nullptr, src.root, condensed_node_mappings);
        
        for (size_t i = 0; i < reads.size(); ++i) {
            std::sort(reads[i]->mutations.begin(), reads[i]->mutations.end());
            remaining_reads.emplace(i);   
        }

        /* map entire set of reads onto tree*/
        this->cartesian_map();
    }
    
    void analyze() {
        for (;!step(TOP_N);) {
            printf("\n");
        }
    }
    
    std::vector<MAT::Node*> get_peaks() {
        return peaks;
    }

    std::vector<MAT::Node*> get_neighbors() {
        return neighbors;
    }
};

//Main peak search algorithm
void analyzeReads(const MAT::Tree &T_ref, const MAT::Tree &T, const std::string &ref_seq, std::unordered_map<size_t, struct read_info*> &read_map, tbb::concurrent_hash_map<MAT::Node*, double> &node_score_map, const std::vector<std::string> &vcf_samples, const std::string &barcode_file, const std::string &read_mutation_depth_vcf, const std::string &condensed_nodes_csv) {
    timer.Start();
    int top_n = 25, prohibited_dist_thresh = 1, neighbor_dist_thresh = 4, neighbor_peaks_thresh = 100, tree_increment, tree_range, tree_overlap = 0, range_factor = 1;
    std::vector<MAT::Node*> peak_nodes, curr_peak_nodes, prohibited_nodes, neighbor_nodes, curr_neighbor_nodes;
    std::vector<size_t> remaining_reads(read_map.size()), remove_reads;
    MAT::Tree T_condensed;
    std::unordered_map<MAT::Node*, std::vector<MAT::Node*>> condensed_node_mappings;
    
    //Decide tree_overlap and tree_range based on range of reads
    for (const auto& rm: read_map) {
        int curr_range = rm.second->end - rm.second->start + 1;
        if (curr_range > tree_overlap)
            tree_overlap = curr_range;
    }
    tree_range = 2 * tree_overlap;
    
    //Create smaller tree with only sites covered by reads 
    createCondensedTree(T.root, read_map, condensed_node_mappings, T_condensed);

    std::cout << "Using fast placement" << std::endl;
    
    std::vector<struct read_info*> read_vector;
    std::transform(read_map.begin(), read_map.end(), std::back_inserter(read_vector), [](const auto &x) { return x.second; });

    AuxManager am{T_condensed, condensed_node_mappings, read_vector};
    am.analyze();
    peak_nodes = am.get_peaks();
    neighbor_nodes = am.get_neighbors();

    std::cout << "Finished fast placement in " << timer.Stop() / 1000 <<  " seconds " << std::endl;
    
    //Verify Recovery of Input Samples
    printf("MUTATION DISTANCE ORIG:\n"); 
    timer.Start();

    std::unordered_set<int> site_read_map;
    for (size_t i = 0; i < read_map.size(); i++) {
        auto rp = read_map.find(i)->second;
        for (int j = rp->start; j <= rp->end; j++)
            site_read_map.insert(j);
    }

    for (auto sample: vcf_samples) {
        int min_dist = 100000;
        auto sample_mutations = getMutations(T_ref, sample);
        //Remove mutations from sample_mutations that are not present in site_read_map
        auto mut_itr = sample_mutations.begin();
        while (mut_itr != sample_mutations.end()) {
            if (site_read_map.find(mut_itr->position) == site_read_map.end())
                mut_itr = sample_mutations.erase(mut_itr);
            else
                mut_itr++;
        }

        MAT::Node* best_node = NULL;
        for (const auto &pn: peak_nodes) {
            for (const auto &node: condensed_node_mappings.find(pn)->second) {
                //Getting node_mutations from the Tree
                auto node_mutations = getMutations(T_ref, node->identifier);
                //Remove mutations from node_mutations that are not present in site_read_map
                auto mut_itr = node_mutations.begin();
                while (mut_itr != node_mutations.end()) {
                    if (site_read_map.find(mut_itr->position) == site_read_map.end())
                        mut_itr = node_mutations.erase(mut_itr);
                    else
                        mut_itr++;
                }

                int curr_dist = mutationDistance(sample_mutations, node_mutations);
                if (curr_dist < min_dist) {
                    min_dist = curr_dist;
                    best_node = pn;
                }
            }
        }
        printf("Node: %s, Closest_node: %s, mutation_distance: %d\n", sample.c_str(), best_node->identifier.c_str(), min_dist);
    }
    fprintf(stderr,"Mutation Distance Verification took %ld sec\n\n", (timer.Stop() / 1000));
    
    //ADD neighbor_nodes to peak_nodes
    peak_nodes.reserve(peak_nodes.size() + neighbor_nodes.size());
    peak_nodes.insert(peak_nodes.end(), neighbor_nodes.begin(), neighbor_nodes.end());
    neighbor_nodes.clear();
    
    generateFilteringData(T, T_condensed, condensed_node_mappings, ref_seq, peak_nodes, read_map, barcode_file, read_mutation_depth_vcf, condensed_nodes_csv);
    condensed_node_mappings.clear();
    MAT::clear_tree(T_condensed);
}