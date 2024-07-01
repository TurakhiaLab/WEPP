#include "wbe.hpp"

constexpr int NUM_RANGE_TREES = 25;
constexpr int MAX_PEAKS = 500;
constexpr int FINAL_PEAKS = 85;
constexpr int TOP_N = 25;
constexpr double READ_DIST_FACTOR_THRESHOLD = 0.5 / 100;
// number of range bins for read distribution
constexpr int NUM_RANGE_BINS = 60;
constexpr int MAX_NEIGHBORS = 150;
constexpr int MAX_PEAK_PEAK_MUTATION = 2;
constexpr int MAX_PEAK_NONPEAK_MUTATION = 4;
constexpr int MAX_CACHED_MULTIPLICITY_SIZE = 1024;
// take scores from top 1 to 0.8 * top 1
constexpr double SINGLE_BATCH_CUTOFF = 0.8;
constexpr double SCORE_EPSILON = 1e-6;
constexpr double ERROR_RATE = 0.02;
constexpr double EM_EPSILON = 0.0;
constexpr double MIN_PROPORTION = 0.5 / 100;

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
    //Read samples.vcf to check how close are peaks to samples 
    std::vector<std::string> vcf_samples;
    readSampleVCF(vcf_samples, vcf_filename_samples);
    //Get the input reads data
    std::unordered_map<size_t, struct read_info*> read_map;
    std::unordered_map<std::string, std::vector<std::string>> reverse_merge;
    load_reads_from_proto(proto_reads, read_map, reverse_merge);

    /*
    std::vector<int> freqs;
    int sum = 0;
    for (size_t i = 0; i < read_map.size(); ++i) {
        freqs.emplace_back(read_map[i]->degree);
        sum += read_map[i]->degree;
    }
    std::sort(freqs.begin(), freqs.end(), std::greater<int>());
    int curr = 0;
    int percent = 0;
    for (size_t i = 0; i < freqs.size(); ++i) {
        curr += freqs[i];
        while (curr * 100 / sum > percent) {
            printf("Percent %d happens at %zu\n", percent, i);
            percent++;
        }
    }
    */
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


    
    //Get haplotype abundances and condensed node names
    timer.Start();
    std::unordered_map<std::string, std::vector<std::string>> condensed_nodeNames_map;
    readCSV(condensed_nodeNames_map, condensed_nodes_csv);
    
    //Get Freyja Lineages
    // readCSV(freyja_lineage_abun_map, freyja_lineage_csv_filename);

    //CREATE condensed tree for neighbor lineage search
    MAT::Tree T_condensed;
    std::unordered_map<MAT::Node*, std::vector<MAT::Node*>> condensed_node_mappings;
    createCondensedTree(T.root, read_map, condensed_node_mappings, T_condensed);
    
    analyzeReads(
        T, 
        ref_seq, 
        read_map, 
        vcf_samples, 
        condensed_nodes_csv,
        barcode_file,
        read_mutation_depth_vcf
    );
    // std::string command = "python src/WBE/peaks_filtering.py " + vm["output-files-prefix"].as<std::string>() + " " + std::string(dir_prefix) + " 100";
    // int result = std::system(command.c_str());
    // if (result)
    //     fprintf(stderr, "\nCannot run peak_filtering.py\n");
}

struct AuxNode {
    // raw score
    double score;
    // 'f' score, where it rewards
    // nodes who's corresponding reads come from a 
    // wide variety of genome places
    double dist_divergence;
    // based on uncompressed tree
    int leaf_count;
    // whether or not it has been selected 
    // as a peak or neighbor
    bool mapped;
    // based on uncompressed tree
    bool is_leaf;

    /* children and mutations  */
    AuxNode* parent;
    /* muts from root to here */
    std::vector<MAT::Mutation> stack_muts;

    // original tree
    std::string id;
    MAT::Node* condensed_source;

    std::vector<AuxNode*> selected_neighbors;

    std::vector<AuxNode*> children;
    // note: does not handle overlapping ranges
    // instead, the i'th bucket corresponds to the interval
    // i * max_genome_size / buckets to (i + 1) * max_genome_size / buckets
    std::array<int, NUM_RANGE_BINS> mapped_read_counts;

    bool has_mutations_in_range(int start, int end) {
        MAT::Mutation search;
        search.position = start;
        auto it = std::lower_bound(condensed_source->mutations.begin(), condensed_source->mutations.end(), search);
        return it != condensed_source->mutations.end() && it->position <= end;
    }

    // given a range and reference mutation list
    // tells the (genome) positions where the this mutation list differs from reference
    // min_pos and max_pos are genome positions
    std::vector<int> mutations(std::vector<MAT::Mutation> const& comp, int min_pos, int max_pos) {
        std::vector<int> muts;

        const int unknown_nuc = 0b1111;

        MAT::Mutation search;
        search.position = min_pos;
        int i = std::lower_bound(stack_muts.begin(), stack_muts.end(), search) - stack_muts.begin();
        search.position = max_pos;
        int last_i = std::upper_bound(stack_muts.begin(), stack_muts.end(), search) - stack_muts.begin();
        int j = 0;

        while (i < last_i || j < (int) comp.size()) {
            if (i == last_i) {                
                if (comp[j].mut_nuc != unknown_nuc) {
                    muts.push_back(comp[j].position);
                }
                ++j;
            }
            else if (stack_muts[i].position < min_pos) {
                ++i;
            }
            else if (stack_muts[i].position > max_pos) {
                assert(0);
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
                if (comp[j].mut_nuc != unknown_nuc) {
                    muts.push_back(comp[j].position);
                }
                ++j;
            }
            else if (stack_muts[i].position == comp[j].position && stack_muts[i].mut_nuc != comp[j].mut_nuc && stack_muts[i].mut_nuc != unknown_nuc) {
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

    int mutation_distance(struct read_info* other) {
        return mutation_distance(other->mutations, other->start, other->end);
    }

    int mutation_distance(AuxNode * other) {
        return mutation_distance(other->stack_muts, 0, INT_MAX);
    }

    double full_score() {
        return score * sqrt(dist_divergence);
    }
};

struct RangedNode {
    AuxNode* root;
    /* the original aux nodes that correspond to this range tree */
    std::vector<AuxNode*> sources;

    /* arena indices of ranged node children */
    std::vector<int> children;

    std::vector<int> mutations(std::vector<MAT::Mutation> const& comp, int min_pos, int max_pos) {
        return root->mutations(comp, min_pos, max_pos);
    }
};

/* used to sort auxnodes by their score */
struct ScoreComparator {
    bool operator() (AuxNode* const& left, AuxNode* const& right) const {
        double const effective_left = left->full_score(); 
        double const effective_right = right->full_score(); 

        if (abs(effective_left - effective_right) > SCORE_EPSILON) {
            return effective_left > effective_right;
        }
        else if (left->leaf_count != right->leaf_count) {
            return left->leaf_count > right->leaf_count;
        }
        else {
            return left->id > right->id;
        }
    }
};

/* used to sort auxnodes by their mutation list */
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
    using mutex_t = tbb::queuing_mutex;

    // list of all the tree nodes
    // whenever we say arena index
    // we mean this list
    std::vector<AuxNode> arena;

    // list of all the range tree nodes
    // whenever we say range tree arena index, it refers
    // to this list
    std::vector<RangedNode> ranged_arena;

    // maps ranges to RangeTree roots (indices)
    // key is a range (say 1000 to 1400)
    // value is the range tree arena index
    std::map<std::pair<int, int>, int> ranged_root_map;
    int num_reads;
    int read_distribution_bin_size;
    std::array<int, NUM_RANGE_BINS> true_read_counts;
    std::array<double, NUM_RANGE_BINS> true_read_distribution;

    const std::vector<read_info *>& reads;
    // given a read index
    // what is its maximum parismony score?
    // i.e max_parismony[i] = best parsimony of ith read
    std::vector<int> max_parismony;
    // given a read index
    // how many epp positions does it map to
    // parsimony_multiplicity[i] = epps of ith read
    std::vector<int> parsimony_multiplicity;
    
    /* if epp_positions_cache[i].size() != multiplicity[i], that means not cached */
    /* since there's too many. Also, values are arena indices */
    // epp_positions_cache[i] = either empty or the arena indices of the epps
    std::vector<std::set<int>> epp_positions_cache;
    
    /* mappable nodes based on their scores */
    /* a node is 'mappable' if has not been marked as a peak or neighbor */
    /* while a fibbonacci heap is more approrpriate, it's overhead is too high*/
    /* so we just sort when necessary */
    std::vector<AuxNode*> current_nodes;

    /* indices of unmapped reads */
    std::set<int> remaining_reads;

    /* map of read_index to paired peak node */
    std::map<int, AuxNode*> correspondence;

    // selected peaks and neighbors
    std::set<AuxNode*, MutationComparator> selected_peaks, selected_neighbors;
    
    // peaks and neighbors in vector form, and using MAT::Node
    std::vector<MAT::Node*> peaks, neighbors;

    /* For a fixed read, we are given the positions where the parent node */
    /* differs from the read, and the positions where the current node */
    /* differs from the read. This function counts how many difference the parent has */
    /* that the child does not. This is NOT the same as the absolute difference in vector sizes */
    /* since in this function, positions the child has that the parent does not */
    /* are of no relevance */
    /* (this is to match reference behavior: can only map if one or mutation positions was removed */
    /* even if total size increases) */
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

    /* maps a single read to entire tree, caching aliveness */
    /* parent_locations is the positions where the parent has different mutations than the read */
    void single_read_tree(RangedNode* curr, std::vector<int> const & parent_locations, struct read_info const*read, std::vector<RangedNode*> &max_nodes, int &max_val) {
        // positions where this node differs from the read */
        std::vector<int> const my_locations = curr->mutations(read->mutations, read->start, read->end);

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

        // if (!curr->root->is_leaf || reductions) {
        if (parsimony < max_val)
        {
            max_val = parsimony;
            max_nodes = {curr};
        }
        else if (parsimony == max_val)
        {
            max_nodes.emplace_back(curr);
        }

        for (int child: curr->children) {
            single_read_tree(&ranged_arena[child], my_locations, read, max_nodes, max_val);
        }
    }

    RangedNode* find_range_tree_for(read_info const * read) {
        std::pair<int, int> search = {read->start, INT_MAX};
        auto it = ranged_root_map.upper_bound(search);
        do {
            it = std::prev(it);
            if (read->start >= it->first.first && read->end <= it->first.second) {
                return &ranged_arena[it->second];
            }
        } while (it != ranged_root_map.begin());

        assert(0);
        return nullptr;
    }

    /* maps a single read to entire tree, caching aliveness */
    /* and finding the epp positions of a given read */
    /* parent_locations is the positions where the parent has different mutations than the read */
    /* max_val corresponds to max parismony */
    /* max indices correspond to the indices of the epp nodes */
    void single_read_tree(struct read_info const* read, std::set<int> &max_indices, int &max_val) {
        std::vector<int> empty_mutation_list;
        std::vector<RangedNode*> max_nodes;
        RangedNode *root = find_range_tree_for(read);
        single_read_tree(root, empty_mutation_list, read, max_nodes, max_val);

        for (RangedNode* rnode: max_nodes) {
            for (AuxNode *src: rnode->sources) {
                max_indices.emplace(src - &arena[0]);
            }
        }
    }

    double node_score(int parsimony, int epps, int degree) {
        return (double) degree / ((1 + parsimony) * epps);
        // return (double) degree / (epps * (1 + parsimony * parsimony));
    }

    /* map all reads to all nodes, maintaining caches if not too large */
    void cartesian_map() {
        mutex_t my_mutex;

        epp_positions_cache.resize(reads.size());
        max_parismony.resize(reads.size());
        parsimony_multiplicity.resize(reads.size());

        tbb::parallel_for(tbb::blocked_range<size_t>(0, reads.size()), 
            [&](tbb::blocked_range<size_t> k) {
            for (size_t r = k.begin(); r != k.end(); ++r) {
                std::set<int> max_indices;
                int max_val = INT32_MAX; /* really want minimum value */

                single_read_tree(reads[r], max_indices, max_val);

                double const delta = node_score(max_val, max_indices.size(), reads[r]->degree);
                {
                    mutex_t::scoped_lock my_lock{my_mutex};
                    for (int arena_index: max_indices) {
                        arena[arena_index].score += delta;
                        // normalization done at the end
                        arena[arena_index].mapped_read_counts[reads[r]->start / this->read_distribution_bin_size] += reads[r]->degree;
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
            /* calculate divergence */
            double divergence = 0;
            for (size_t j = 0; j < NUM_RANGE_BINS; ++j) {
                double const proportion = (double) arena[i].mapped_read_counts[j] / true_read_counts[j];
                if (proportion > READ_DIST_FACTOR_THRESHOLD) {
                    divergence += 1.0 / NUM_RANGE_BINS;
                }
            }
            arena[i].dist_divergence = divergence;

            this->current_nodes.emplace_back(&arena[i]);
        }

        /* initial sort into scores */
        std::sort(current_nodes.begin(), current_nodes.end(), ScoreComparator());
    }

    /* find all current reads that map to this node */
    std::vector<int> find_correspondents(AuxNode *node) {
        int const arena_index = node - &arena[0];
        std::vector<int> correspondents;

        mutex_t my_mutex;

        std::vector<int> remaining(remaining_reads.begin(), remaining_reads.end());
        tbb::parallel_for(tbb::blocked_range<size_t>(0, remaining.size()),
            [&](tbb::blocked_range<size_t> r) {
                for (size_t i = r.begin(); i != r.end(); ++i) {
                    size_t const read = remaining[i]; 
                    /* if cached correspond, use that */
                    if (epp_positions_cache[read].find(arena_index) != epp_positions_cache[read].end()) {
                        mutex_t::scoped_lock lock{my_mutex};
                        correspondents.push_back(read); 
                    }
                    else if (epp_positions_cache[read].size() == (size_t) parsimony_multiplicity[read]) {
                        /* cached and not found */
                        continue;
                    }
                    else {
                        /* recompute to see if correspondent */
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
                            mutex_t::scoped_lock lock{my_mutex};
                            correspondents.push_back(read);
                        }
                    }
                } 
            }
        );
        std::sort(correspondents.begin(), correspondents.end());

        return correspondents;
    }

    /* removes a singular read's contribution from current tree */
    void remove_reads_effects(int index, mutex_t * mutex) {
        std::set<int> recalcuated;
        std::set<int> *epps;
        /* if cached, just take that */
        if (this->epp_positions_cache[index].size() == (size_t) parsimony_multiplicity[index]) {
            epps = &this->epp_positions_cache[index];
        }
        else {
            /* recalculate on (almost) entire tree, generally the most optimal */
            int max_val = INT32_MAX;
            single_read_tree(reads[index], recalcuated, max_val);
            epps = &recalcuated;
        }

        /* use ORIGINAL size */
        mutex_t::scoped_lock lock{*mutex};
        double const delta = node_score(max_parismony[index], parsimony_multiplicity[index], reads[index]->degree);
        for (int arena_index: (*epps)) {
            AuxNode *node = &arena[arena_index];
            if (node->mapped) {
                continue;
            }
            
            node->score -= delta;
        }
    }

    void singular_step(AuxNode* node) {
        /* 1. map remaining reads onto this node (finding correspondents) */
        std::vector<int> correspondents = find_correspondents(node);

        /* 2. map all said reads onto entire remaining set, adjusting deltas */
        mutex_t my_mutex;

        tbb::parallel_for(tbb::blocked_range<size_t>(0, correspondents.size()),
            [&](tbb::blocked_range<size_t> range) {
                for (size_t i = range.begin(); i < range.end(); ++i) {
                    this->remove_reads_effects(correspondents[i], &my_mutex);
                }
            }
        );

        /* 3. mark correspondence */ 
        for (int correspondent: correspondents) {
            correspondence[correspondent] = node;
            this->remaining_reads.erase(correspondent);
        }
    }

    /* can both aux nodes be sent in the same batch? */
    bool valid_two_tops(AuxNode *a, AuxNode *b) {
         return a->mutation_distance(b) > MAX_PEAK_PEAK_MUTATION;
    }

    /* adds possible neighbors of subtree to the given set if the mutation distance is less than a threshold */
    void dfs_possible_neighbors(AuxNode* pivot, AuxNode *curr, std::set<AuxNode*, ScoreComparator>& s, int radius) {
        if (pivot->mutation_distance(curr) > radius) {
            return;
        }

        if (!curr->mapped) {
            s.insert(curr);
        }

        for (AuxNode* child: curr->children) {
            dfs_possible_neighbors(pivot, child, s, radius);
        }
    }

    /* find all possible neighbors within a radius from a given node */
    void possible_neighbors(AuxNode* pivot, std::set<AuxNode*, ScoreComparator>& s, int radius) {
        AuxNode *curr = pivot;
        while (curr->parent && pivot->mutation_distance(curr->parent) <= radius) {
            curr = curr->parent;
        }

        dfs_possible_neighbors(pivot, curr, s, radius);
    }

    void mark_as_peaks(std::vector<AuxNode*> const& nodes) {
        for (AuxNode* const& n: nodes) {
            selected_peaks.insert(n);
            peaks.emplace_back(n->condensed_source);
        }
    }

    /* clears neighbors (and self) from map */
    void clear_neighbors(std::vector<AuxNode*> const &nodes) {
        /* 1. mark as peak */
        this->mark_as_peaks(nodes);

        std::vector<AuxNode*> added;
        /* 2. find candidates within radius (taking only top n for each pivot) */
        /* buffer the selected neighbors into added */
        for (AuxNode *pivot: nodes) {
            std::set<AuxNode*, ScoreComparator> multisource_radius;
            possible_neighbors(pivot, multisource_radius, MAX_PEAK_NONPEAK_MUTATION);

            pivot->selected_neighbors.clear();
            int i = 0;
            for (AuxNode * node : multisource_radius)
            {
                bool const found = selected_peaks.find(node) != selected_peaks.end() || selected_neighbors.find(node) != selected_neighbors.end();
                if (!found) {
                    node->mapped = true;
                    added.emplace_back(node);
                    neighbors.emplace_back(node->condensed_source);
                    pivot->selected_neighbors.emplace_back(node);
                    if (++i == MAX_NEIGHBORS)
                    {
                        break;
                    }
                }
            }
        }

        /* 3. for all nodes within small radius, forcefuly remove even if more than 100 */        
        for (AuxNode *pivot: nodes) {
            std::set<AuxNode*, ScoreComparator> multisource_radius;
            possible_neighbors(pivot, multisource_radius, MAX_PEAK_PEAK_MUTATION);

            for (AuxNode * node : multisource_radius)
            {
                node->mapped = true;
            }
        }

        selected_neighbors.insert(added.begin(), added.end());
    }

    /* matches currently top_n nodes with their reads */
    /* returns if finished or not */
    bool step() {
        assert(!current_nodes.empty() && !remaining_reads.empty());

        /* to get exactly same behavior as before, have a single considering set */
        std::vector<AuxNode*> consideration;

        /* get top_n (under some special constraints) */ 
        auto it = current_nodes.begin();
        double const min_score = (*it)->full_score();

        /* no available peaks */
        if (min_score < SCORE_EPSILON) {
            return true;
        }

        /* while: */
        /* 1. we have current_nodes left */
        /* 2. the score is all the samae */
        /* 3. we have not selected 25 yet */
        /* 4. we have not exceeded max peaks */
        for (int i = 0; 
                it != current_nodes.end() && 
                abs((*it)->full_score() - min_score) < SCORE_EPSILON &&
                consideration.size() < (size_t) TOP_N &&
                consideration.size() + selected_peaks.size() < MAX_PEAKS;
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
                printf("%.9f raw %.9f divergence id: %s\n", (*it)->score, (*it)->dist_divergence, (*it)->id.c_str());
            }
        }

        /* clear neighbors of selected peaks (marking them as mapped as well) */
        this->clear_neighbors(consideration);

        for (AuxNode *&node: consideration) {
            this->singular_step(node);
        }

        /* resort current_nodes (removing mapped nodes as well) */
        current_nodes.erase(
            std::remove_if(current_nodes.begin(), current_nodes.end(),
                [](AuxNode* const& curr) {
                    return curr->mapped || curr->score <= SCORE_EPSILON;
                }
            ),
            current_nodes.end()
        );
        std::sort(current_nodes.begin(), current_nodes.end(), ScoreComparator());

        return selected_peaks.size() >= MAX_PEAKS || remaining_reads.empty() || current_nodes.empty();
    }

    /* construct auxiliary tree from mat_tree */
    AuxNode* from_mat_tree(AuxNode* parent, MAT::Node* node, std::unordered_map<MAT::Node*, std::vector<MAT::Node*>>& condensed_node_mappings) {
        this->arena.emplace_back();
        AuxNode *ret = &this->arena.back();
        ret->parent = parent;
        ret->mapped = false;
        ret->is_leaf = condensed_node_mappings.at(node).front()->is_leaf(); 
        ret->score = 0;
        ret->stack_muts = {};
        ret->id = node->identifier;
        ret->condensed_source = node;
        ret->leaf_count = getNumLeaves(condensed_node_mappings, node);
        /* flatten mutation list */
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

    int build_range_tree(int parent, AuxNode* curr, int start, int end) {
        int ret;
        if (parent == -1 || curr->has_mutations_in_range(start, end)) {
            ret = ranged_arena.size();
            ranged_arena.emplace_back();
            ranged_arena[ret].root = curr;
            ranged_arena[ret].sources = {curr};
        }
        else {
            ret = parent;
            ranged_arena[parent].sources.push_back(curr);
        }

        for (AuxNode* child: curr->children) {
            int const sub = build_range_tree(ret, child, start, end);
            if (sub != ret) {
                ranged_arena[ret].children.emplace_back(sub);
            }
        }

        return ret;
    }

    void build_range_trees() {
        // get ranges (heuristic based approach)
        using range = std::pair<int, int>;
        std::vector<range> read_ranges;
        std::transform(reads.begin(), reads.end(), std::back_inserter(read_ranges),
            [](read_info * const &read) {
                return std::make_pair(read->start, read->end);
            }
        );
        std::sort(read_ranges.begin(), read_ranges.end());

        int const num = std::min((int) read_ranges.size(), NUM_RANGE_TREES);
        for (int i = 0; i < num; ++i) {
            int s = read_ranges.size() * i / num;
            int e = read_ranges.size() * (i + 1) / num;
            int read_start = read_ranges[s].first;
            int read_end = 0;
            for (int j = s; j < e; ++j) {
                read_end = std::max(read_end, read_ranges[j].second);
            }

            // there is technically a chance we can have a collision 
            // where multiple range trees are for the exact same range
            // (to see this, if all reads are the exact same range
            // then all range trees will be the exact same)
            // However, this is not an issue even in the case that there are duplicates, 
            // it only reduces performance, it doesn't actually make anything break
            range r{read_start, read_end};
            if (this->ranged_root_map.find(r) == ranged_root_map.end()) {
                this->ranged_root_map[r] = build_range_tree(-1, &arena[0], read_start, read_end);
            }
        }

        this->read_distribution_bin_size = (read_ranges.back().first + NUM_RANGE_BINS - 1) / NUM_RANGE_BINS;
        this->num_reads = 0;
        std::fill(true_read_counts.begin(), true_read_counts.end(), 0);
        for (const auto& r: reads) {
            true_read_counts[r->start / read_distribution_bin_size] += r->degree;
            this->num_reads += r->degree;
        }
        for (size_t i = 0; i < NUM_RANGE_BINS; ++i) {
            true_read_distribution[i] = (double) true_read_counts[i] / this->num_reads;
        }
    }

    void em_postprocess(bool include_neighbors) {
        current_nodes.clear();
        tbb::parallel_for(tbb::blocked_range<size_t>(0, arena.size()),
            [&](tbb::blocked_range<size_t> range) {
                for (size_t i = range.begin(); i < range.end(); ++i) {
                    arena[i].mapped = false;
                    arena[i].score = 0;
                    arena[i].dist_divergence = 1;
                }
            } 
        );

        std::vector<AuxNode*> subset;
        subset.insert(subset.end(), selected_peaks.begin(), selected_peaks.end());
        if (include_neighbors) {
            for (AuxNode* peak: selected_peaks) {
                subset.insert(subset.end(), peak->selected_neighbors.begin(), peak->selected_neighbors.end());
            }
        }
        current_nodes.insert(current_nodes.end(), subset.begin(), subset.end());

        std::vector<std::vector<double>> q(subset.size(), std::vector<double>(reads.size()));
        tbb::parallel_for(tbb::blocked_range<size_t>(0, subset.size()),
             [&](tbb::blocked_range<size_t> range) {
                for (size_t i = range.begin(); i < range.end(); ++i) {
                    for (size_t j = 0; j < reads.size(); ++j) {
                        double dist = subset[i]->mutation_distance(reads[j]);
                        q[i][j] = std::pow(ERROR_RATE, dist) * std::pow(1 - ERROR_RATE, reads[j]->end - reads[j]->start + 1 - dist);
                    }
                }
             }
        );

        double prev = 0, curr = 0;
        int max_it = 50;
        std::vector<double> p(subset.size(), 1.0 / subset.size());
        do {
            prev = curr;

            std::vector<double> denom(reads.size());
            tbb::parallel_for(tbb::blocked_range<size_t>(0, reads.size()), 
                [&](tbb::blocked_range<size_t> range) {
                    for (size_t j = range.begin(); j < range.end(); ++j) {
                        double sum = 0;
                        for (size_t l = 0; l < subset.size(); ++l) {
                            sum += p[l] * q[l][j];
                        }
                        denom[j] = sum;
                    }
                }
            );

            tbb::parallel_for(tbb::blocked_range<size_t>(0, subset.size()), 
                [&](tbb::blocked_range<size_t> range) {
                    for (size_t i = range.begin(); i < range.end(); ++i) {
                        double sum = 0;
                        size_t count = 0;
                        for (size_t j = 0; j < reads.size(); ++j) {
                            sum += reads[j]->degree * p[i] * q[i][j] / denom[j];
                            count += reads[j]->degree;
                        }
                        // don't even need multiple p
                        p[i] = sum / count;
                    }
                } 
            );

            std::vector<double> logl(reads.size());
            tbb::parallel_for(tbb::blocked_range<size_t>(0, reads.size()), 
                [&](tbb::blocked_range<size_t> range) {
                    for (size_t j = range.begin(); j < range.end(); ++j) {
                        double sum = 0;
                        for (size_t l = 0; l < subset.size(); ++l) {
                            sum += p[l] * q[l][j];
                        }
                        logl[j] = std::log(sum);
                    }
                } 
            );

            curr = std::accumulate(logl.begin(), logl.end(), 0.0);
        } while (abs(curr - prev) > EM_EPSILON && --max_it);

        for (size_t i = 0; i < p.size(); ++i) {
            subset[i]->score = p[i];
        }
        std::sort(subset.begin(), subset.end(), ScoreComparator());
        auto it = subset.begin();
        double total = 0;
        while (it != subset.end() && (*it)->score / ((*it)->score + total) >= MIN_PROPORTION) {
            total += (*it)->score;
            ++it;
        }
        subset.erase(it, subset.end());

        this->peaks = {};
        this->selected_peaks = this->selected_neighbors = {};
        this->mark_as_peaks(subset);

        if (include_neighbors) {
            std::cout << "****Secondary EM (" << subset.size() << ")****" << std::endl;
        }
        else {
            std::cout << "****Initial EM (" << subset.size() << ")****" << std::endl;
        }

        for (size_t i = 0; i < subset.size(); ++i) {
            std::cout << "EM Selected " << subset[i]->id << " with abundance " << std::fixed << std::setprecision(9) << subset[i]->score / total << std::endl;
        }
    }

    void postprocess() {
        double const sigma = 4.0;
        double const coeff = log(1 / (sigma * sqrt(2 * M_PI)));

        current_nodes.clear();
        tbb::parallel_for(tbb::blocked_range<size_t>(0, arena.size()),
            [&](tbb::blocked_range<size_t> range) {
                for (size_t i = range.begin(); i < range.end(); ++i) {
                    arena[i].mapped = false;
                    arena[i].score = 0;
                    arena[i].dist_divergence = 1;
                }
            } 
        );

        std::vector<AuxNode*> subset;
        subset.insert(subset.end(), selected_peaks.begin(), selected_peaks.end());
        for (AuxNode* node: selected_peaks) {
            subset.insert(subset.end(), node->selected_neighbors.begin(), node->selected_neighbors.end());
        }
        current_nodes.insert(current_nodes.end(), subset.begin(), subset.end());
        
        std::vector<std::vector<double>> parsimony(reads.size(), std::vector<double>(subset.size()));
        tbb::parallel_for(tbb::blocked_range<size_t>(0, subset.size()),
             [&](tbb::blocked_range<size_t> range) {
                for (size_t i = range.begin(); i < range.end(); ++i) {
                    for (size_t j = 0; j < reads.size(); ++j) {
                        double dist = subset[i]->mutation_distance(reads[j]) / sigma;
                        parsimony[j][i] = reads[j]->degree * (coeff - dist * dist);
                        subset[i]->score += parsimony[j][i];
                    }
                }
             }
        );

        selected_peaks.clear();
        selected_neighbors.clear();
        peaks.clear();
        neighbors.clear();

        std::vector<AuxNode*> consideration;

        for (size_t q = 0; q < subset.size() && q < FINAL_PEAKS; ++q) {
            std::sort(current_nodes.begin(), current_nodes.end(), ScoreComparator());

            current_nodes[0]->mapped = true;
            consideration.emplace_back(current_nodes[0]);

            tbb::parallel_for(tbb::blocked_range<size_t>(0, subset.size()),
                [&](tbb::blocked_range<size_t> range) {
                    for (size_t i = range.begin(); i < range.end(); ++i) {
                        subset[i]->score = 0;
                    }
                }
            );

            std::vector<double> comp(reads.size());
            tbb::parallel_for(tbb::blocked_range<size_t>(0, reads.size()),
                [&](tbb::blocked_range<size_t> range) {
                    for (size_t i = range.begin(); i < range.end(); ++i) {
                        double dist = current_nodes[0]->mutation_distance(reads[i]) / sigma;
                        comp[i] = reads[i]->degree * (coeff - dist * dist);
                    }
                }
            );

            tbb::parallel_for(tbb::blocked_range<size_t>(0, subset.size()),
                [&](tbb::blocked_range<size_t> range) {
                    for (size_t i = range.begin(); i < range.end(); ++i) {
                        for (size_t j = 0; j < reads.size(); ++j) {
                            parsimony[j][i] = std::max(parsimony[j][i], comp[j]);
                            subset[i]->score += parsimony[j][i];
                        }
                    }
                }
            );

            current_nodes.erase(
                std::remove_if(
                    current_nodes.begin(),
                    current_nodes.end(),
                    [](AuxNode* const &node) -> bool {
                        return node->mapped;
                    }
                ),
                current_nodes.end()
            );
        }

        /* clear neighbors */
        this->mark_as_peaks(consideration);
    }

    void kmeans_postprocess(double min_abundance) {
        constexpr int max_it = 15, explore_rad = 5;
        constexpr double stay_factor = 0.75;
        std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());
        
        const size_t n = (size_t) (1 + 1 / min_abundance);
        
        tbb::parallel_for(tbb::blocked_range<size_t>(0, arena.size()),
            [&](tbb::blocked_range<size_t> range) {
                for (size_t i = range.begin(); i < range.end(); ++i) {
                    arena[i].mapped = false;
                    arena[i].score = 0;
                    arena[i].dist_divergence = 1;
                }
            } 
        );

        // initially take results from WEPP
        std::vector<AuxNode*> selected(selected_peaks.begin(), selected_peaks.end());
        // std::vector<AuxNode*> selected;
        // for (size_t j = 0; j < n; ++j) {
            // selected.push_back(&arena[rng() % arena.size()]);
        // }
        std::cout << " Selected Size " << selected.size() << std::endl;

        size_t total_reads_degree = 0;
        for (size_t i = 0; i < reads.size(); ++i) {
            total_reads_degree += reads[i]->degree;
        }

        for (int k = 0; k < max_it; ++k) {
            my_mutex_t mutex;
            std::vector<std::vector<int>> corresponding_reads(selected.size());
            std::vector<std::pair<double, int>> index_set(selected.size());
            for (size_t i = 0; i < selected.size(); ++i) {
                index_set[i].second = i;
            }

            // map all reads to current peaks
            // and create epp sets
            tbb::parallel_for(tbb::blocked_range<size_t>(0, reads.size()), 
                [&](tbb::blocked_range<size_t> range) {
                    for (size_t i = range.begin(); i < range.end(); ++i) {
                        int min_dist = INT_MAX;
                        std::vector<size_t> best;
                        for (size_t j = 0; j < selected.size(); ++j) {
                            int dist = selected[j]->mutation_distance(this->reads[i]);
                            if (dist < min_dist) {
                                min_dist = dist;
                                best = {j};
                            }
                            else if (dist == min_dist) {
                                best.emplace_back(j);    
                            }
                        }

                        // avoiding data race on both rng and corresponding_reads
                        my_mutex_t::scoped_lock _lock{mutex};
                        int ind = best[rng() % best.size()];
                        corresponding_reads[ind].emplace_back(i);
                        index_set[ind].first += (double) reads[ind]->degree / total_reads_degree;
                    }
                }
            );

            std::sort(index_set.begin(), index_set.end());
            // index of previous selection
            std::vector<size_t> seeds;
            size_t j = 0;
            double running_abundance = 0;
            // we delete nodes with too little abundance
            // and split nodes with too high abudance
            // idea is basically that
            // < abundance / 2 = delete
            // > abundance * 2 = split (theoretically to same node)
            // while (j < index_set.size()) {
            //     if (index_set[j].first + running_abundance + 1e-6 >= min_abundance) {
            //         seeds.emplace_back(index_set[j].second);
            //         index_set[j].first -= min_abundance - running_abundance;
            //         running_abundance = 0;
            //     }
            //     else {
            //         running_abundance += index_set[j].first;
            //         index_set[j].first = 0;
            //         ++j;
            //     } 
            // }
            std::vector<AuxNode*> new_selection;

            // debug (ignore running abundance)
            //seeds = std::vector<size_t>(selected.size());
            //std::iota(seeds.begin(), seeds.end(), 0);
            seeds = {};
            for (size_t j = 0; j < index_set.size(); ++j) {
                if (index_set[j].first < min_abundance * 2) {
                    seeds.push_back(index_set[j].second);
                }
                else {
                    seeds.push_back(index_set[j].second);
                    // seeds.push_back(index_set[j].second);
                }
            }
            
            // map n -> n'
            for (size_t i = 0; i < seeds.size(); ++i) {
                std::set<AuxNode *, ScoreComparator> all_neighbors;
                possible_neighbors(selected[seeds[i]], all_neighbors, explore_rad);

                // randomly select weighted neighbor based on scores
                my_mutex_t mutex;

                std::vector<AuxNode*> flat(all_neighbors.begin(), all_neighbors.end());
                std::vector<size_t> distances;
                size_t best_distance = SIZE_MAX;
                for (AuxNode* node: flat) {
                    size_t total_sum = 0;
                    tbb::parallel_for(tbb::blocked_range<size_t>(0, corresponding_reads[seeds[i]].size()), [&](tbb::blocked_range<size_t> range) {
                        size_t sum = 0;
                        for (size_t j = range.begin(); j < range.end(); ++j)
                        {
                            struct read_info* read = this->reads[corresponding_reads[seeds[i]][j]];
                            sum += (size_t) read->degree * node->mutation_distance(read);
                        }

                        my_mutex_t::scoped_lock _lock{mutex};
                        total_sum += sum; 
                    });

                    distances.push_back(total_sum);
                    best_distance = std::min(best_distance, total_sum);
                }

                // lower average -> higher chance of selection
                std::vector<double> prob;
                double normalization = 0;
                for (size_t j = 0; j < flat.size(); ++j) {
                    normalization += exp(-stay_factor * (double) (distances[j] - best_distance) / flat.size());
                    prob.push_back(normalization);
                }
                std::uniform_real_distribution<double> unif(0, normalization);
                double rand = unif(rng);
                double prev = 0;
                AuxNode* best_node = flat.back();
                for (size_t j = 0; j < prob.size(); ++j) {
                    if (prev <= rand && rand < prob[j]) {
                        best_node = flat[j];
                        break;
                    }
                    prev = prob[j];
                }
                new_selection.emplace_back(best_node);
            }

            selected = new_selection;
            std::cout << "Selection Size raw " << selected.size() <<  " unique " << std::set<AuxNode*>(selected.begin(), selected.end()).size() << std::endl;
        }

        selected_peaks = {};
        peaks = {};
        mark_as_peaks(selected);
    }

    void init_from_peaks( std::unordered_map<MAT::Node*, std::vector<MAT::Node*>>& condensed_node_mappings) {
        std::set<std::string> peak_names = {
            "hCoV-19/England/MILK-9E05B3/2020|EPI_ISL_601443|2020-09-20",
            "hCoV-19/Japan/IC-0564/2021|EPI_ISL_792683|2021-01-02",
            "hCoV-19/USA/NY-MSHSPSP-PV24650/2020|EPI_ISL_1300881|2020-12-06",
            "hCoV-19/India/MH-NCCS-P1162000182735/2021|EPI_ISL_1544014|2021-02-27",
            "hCoV-19/Hong",
            "hCoV-19/Australia/QLD2568/2021|EPI_ISL_7190366|2021-12-01",
            "hCoV-19/England/MILK-2DF642C/2021|EPI_ISL_7718520|2021-12-09",
            "hCoV-19/Denmark/DCGC-493190/2022|EPI_ISL_12248637|2022-04-17"
        };
        std::vector<AuxNode*> peaks;
        for (size_t i = 0; i < arena.size(); ++i) {
            for (auto node: condensed_node_mappings[arena[i].condensed_source]) {
                if (peak_names.find(node->identifier) != peak_names.end()) {
                    peaks.emplace_back(&arena[i]);
                }
            }
        }

        clear_neighbors(peaks);
    }

    void init_from_lineages(std::unordered_map<MAT::Node*, std::vector<MAT::Node*>>& condensed_node_mappings)
    {
        std::set<std::string> lineage_names = {
            "BG.5","AY.30","BA.2","Q.1","B.1.526","BA.1.7","XAH","P.1.4","BA.2.23","BA.1","BA.1.17.2","B.1.1.529","BA.1.9","BA.2.33","B.1.617.2","BA.2.15","B.1.1.7","P.1.7","AY.50","BA.1.1.10","BA.2.27","BA.1.5","BA.2.9.4","BA.2.29","P.1.11","AY.17"
        };
        std::vector<AuxNode*> peaks;
        for (size_t i = 0; i < arena.size(); ++i) {
            for (auto node: condensed_node_mappings[arena[i].condensed_source]) {
                auto lineage = node->clade_annotations[1];
                if (lineage_names.find(lineage) != lineage_names.end()) {
                    lineage_names.erase(lineage);
                    peaks.emplace_back(&arena[i]);
                }
            }
        }

        std::cout << " Peaks Size " << peaks.size() << std::endl;
        mark_as_peaks(peaks);
    }

    void init_from_indices()
    {
        std::vector<int> indices{3173879 ,2558943 ,1742551 ,1388539 ,2566532 ,3221605 ,3469427 ,2216187 ,1565422 ,1720997 ,1569271 ,3479779 ,3179800 ,3236899 ,1550036 ,1551669 ,2557539 ,3479749 ,3254747 ,3471804 ,1565423 ,1567175 ,3471805 ,1388478 ,3223973 ,2729994 ,1564004 ,3325633 ,1568871 ,1388594 ,2435135 ,1564054 ,3318665 ,1564293 ,1388524};
        std::vector<AuxNode*> peaks;
        for (size_t i = 0; i < indices.size(); ++i) {
            peaks.emplace_back(&arena[indices[i]]);
        };

        mark_as_peaks(peaks);
    }

    // for freyja
    void dump_barcode(std::vector<AuxNode*> inputs)
    {
        std::set<std::string> mutations;
        std::vector<std::set<std::string>> node_muts;
        for (AuxNode* n: inputs) {
            std::set<std::string> my_muts;
            for (MAT::Mutation mut: n->stack_muts) {
                std::string build = mut.get_string();
                mutations.insert(build);
                my_muts.insert(std::move(build));
            } 
            node_muts.push_back(std::move(my_muts));
        }
        std::vector<std::string> mutation_vec(mutations.begin(), mutations.end());

        std::ofstream outfile("Freyja/data/usher_barcodes.csv");
        for (const std::string& mut: mutations) {
            outfile << "," << mut;
        }
        outfile << std::endl;
        for (size_t i = 0; i < inputs.size(); ++i) {
            AuxNode* n = inputs[i];
            outfile << 'N' << (n - &arena[0]);
            for (const std::string &mut : mutations)
            {
                bool contains = node_muts[i].find(mut) != node_muts[i].end();
                outfile << "," << (contains ? '1' : '0');
            }
            outfile << std::endl;
        }
    }

    void freyja_postprocess() {
        std::vector<AuxNode*> all_selected(this->selected_peaks.begin(), this->selected_peaks.end());
        all_selected.insert(all_selected.end(), this->selected_neighbors.begin(), this->selected_neighbors.end());
        dump_barcode(all_selected);

        if (std::system("bash -c \"source ~/miniconda3/etc/profile.d/conda.sh && conda activate freyja-env && freyja demix Freyja/cwap_variants.tsv Freyja/cwap_depth.tsv --barcodes Freyja/data/usher_barcodes.csv --output Freyja/my_output_latest.txt\"") != 0) {
            std::cerr << "Failed to run freyja" << std::endl;
            return;
        }

        std::vector<AuxNode*> freyja_nodes;

        std::ifstream fin("Freyja/my_output_latest.txt");
        std::string tmp;
        std::getline(fin, tmp);
        std::getline(fin, tmp);
        fin >> tmp;
        // all of the selected ids
        std::getline(fin, tmp);
        std::stringstream ss{tmp};
        std::string index;
        while (ss >> index) {
            // skip past the 'N'
            int ind = atoi(index.c_str() + 1);
            freyja_nodes.push_back(&arena[ind]);
        }


        this->selected_peaks = {};
        this->selected_neighbors = {};
        this->peaks = {};
        this->neighbors = {};
        this->mark_as_peaks(freyja_nodes);
    }

public:
    /* preprocessing */
    AuxManager(const MAT::Tree &src, std::unordered_map<MAT::Node*, std::vector<MAT::Node*>>& condensed_node_mappings, std::vector<read_info*> const& reads) : reads{reads} {
        /* create auxiliary tree, ensure no reallocations */
        this->arena.reserve(mat_tree_size(src.root));
        this->from_mat_tree(nullptr, src.root, condensed_node_mappings);
        this->build_range_trees();
        
        for (size_t i = 0; i < reads.size(); ++i) {
            std::sort(reads[i]->mutations.begin(), reads[i]->mutations.end());
            remaining_reads.emplace(i);   
        }

        /* map entire set of reads onto tree */
        this->cartesian_map();
        // this->init_from_peaks(condensed_node_mappings);
        // this->init_from_lineages(condensed_node_mappings);

        // this->init_from_peaks(condensed_node_mappings);
        // std::vector<AuxNode*> all_nodes(selected_peaks.begin(), selected_peaks.end());
        // all_nodes.insert(all_nodes.end(), selected_neighbors.begin(), selected_neighbors.end());
        // this->dump_barcode(all_nodes);
    }
    
    void analyze() {
        for (;!step();) {
            printf("\n");
        }

        int it = 1; 
        for (int i = 0; i < it; ++i) {
            this->freyja_postprocess();

            if (i < it - 1) {
                std::vector<AuxNode*> selected(this->selected_peaks.begin(), this->selected_peaks.end());
                this->clear_neighbors(selected);
            }
        }
    }
    
    std::vector<MAT::Node*> get_peaks(void) {
        return peaks;
    }

    std::vector<MAT::Node*> get_neighbors(void) {
        return neighbors;
    }
};

//Main peak search algorithm
void analyzeReads(
    const MAT::Tree &T, 
    const std::string &ref_seq, 
    std::unordered_map<size_t, struct read_info*> &read_map,
    const std::vector<std::string> &vcf_samples,
    const std::string &condensed_nodes_csv,
    const std::string &barcode_file,
    const std::string &read_mutation_depth_vcf
) {
    timer.Start();
    std::vector<MAT::Node*> peak_nodes, curr_peak_nodes, prohibited_nodes, neighbor_nodes, curr_neighbor_nodes;
    std::vector<size_t> remaining_reads(read_map.size()), remove_reads;
    MAT::Tree T_condensed;
    std::unordered_map<MAT::Node*, std::vector<MAT::Node*>> condensed_node_mappings;
    
    //Create smaller tree with only sites covered by reads 
    createCondensedTree(T.root, read_map, condensed_node_mappings, T_condensed);

    std::cout << "Using fast placement" << std::endl;
    
    std::vector<struct read_info*> read_vector;
    std::transform(read_map.begin(), read_map.end(), std::back_inserter(read_vector), [](const auto &x) { return x.second; });

    AuxManager am{T_condensed, condensed_node_mappings, read_vector};
    std::cout << "Finished cartesian mapping in " << timer.Stop() / 1000 <<  " seconds " << std::endl;
    am.analyze();
    peak_nodes = am.get_peaks();
    neighbor_nodes = am.get_neighbors();

    std::cout << "Finished fast placement in " << timer.Stop() / 1000 <<  " seconds " << std::endl;
    std::cout << " Selected " << (peak_nodes.size() + neighbor_nodes.size()) << std::endl;
    
    //Verify Recovery of Input Samples
    printf("MUTATION DISTANCE ORIG:\n"); 
    timer.Start();

    // peak_nodes.insert(peak_nodes.end(), neighbor_nodes.begin(), neighbor_nodes.end());

    std::unordered_set<int> site_read_map;
    for (size_t i = 0; i < read_map.size(); i++) {
        auto rp = read_map.find(i)->second;
        for (int j = rp->start; j <= rp->end; j++)
            site_read_map.insert(j);
    }

    double average_distance = 0;
    for (auto sample: vcf_samples) {
        int min_dist = 100000;
        auto sample_mutations = getMutations(T, sample);
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
                auto node_mutations = getMutations(T, node->identifier);
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
        average_distance += (float) min_dist / vcf_samples.size();
        printf("Node: %s, Closest_node: %s, mutation_distance: %d\n", sample.c_str(), best_node->identifier.c_str(), min_dist);
    }
    fprintf(stderr,"Mutation Distance Verification took %ld sec\n\n", (timer.Stop() / 1000));
    printf("Average Distance %f\n", average_distance);
    
    // //ADD neighbor_nodes to peak_nodes
    peak_nodes.reserve(peak_nodes.size() + neighbor_nodes.size());
    peak_nodes.insert(peak_nodes.end(), neighbor_nodes.begin(), neighbor_nodes.end());
    neighbor_nodes.clear();

    // generateFilteringData(T, T_condensed, condensed_node_mappings, ref_seq, peak_nodes, read_map, barcode_file, read_mutation_depth_vcf, condensed_nodes_csv);
}