#include "wbe.hpp"

constexpr int NUM_RANGE_TREES = 25;
constexpr int MAX_PEAKS = 7500;
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
constexpr double EM_EPSILON = 0.1;
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
    // //RUN peaks_filtering
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
        /* 2. the score is relatively the same */
        /* 3. we have not exceeded max peaks */
        for (int i = 0; 
                it != current_nodes.end() && 
                abs((*it)->full_score() / min_score) > SINGLE_BATCH_CUTOFF &&
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

        /* remove associated reads of ONLY the top */ 
        if (!consideration.empty()) {
            this->singular_step(consideration[0]);
        }
        // for (AuxNode *&node: consideration) {
        //     this->singular_step(node);
        // }

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
        int max_it = 15;
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

        this->selected_peaks = this->selected_neighbors = {};
        this->mark_as_peaks(subset);

        if (include_neighbors) {
            std::cout << "****Final EM (" << subset.size() << ")****" << std::endl;
        }
        else {
            std::cout << "****Initial EM (" << subset.size() << ")****" << std::endl;
        }

        for (size_t i = 0; i < subset.size(); ++i) {
            std::cout << "EM Selected " << subset[i]->id << " with abundance " << std::fixed << std::setprecision(9) << subset[i]->score / total << std::endl;
        }
    }

    void postprocess() {
        double const sigma = 3.0;
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
        subset.insert(subset.end(), selected_neighbors.begin(), selected_neighbors.end());
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

        for (size_t q = 0; q < subset.size() && q < MAX_PEAKS; ++q) {
            std::sort(current_nodes.begin(), current_nodes.end(), ScoreComparator());

            current_nodes[0]->mapped = true;
            consideration.emplace_back(current_nodes[0]);

            std::cout << " Working score " << current_nodes[0]->score << " name " << current_nodes[0]->id << std::endl;

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
        this->clear_neighbors(consideration);
    }

    void init_from_peaks( std::unordered_map<MAT::Node*, std::vector<MAT::Node*>>& condensed_node_mappings) {
        std::set<std::string> peak_names = {
            "USA/UT-UPHL-211207411157/2021|OL850427.1|2021-10-23","England/DHSC-CYNGP8M/2022|2022-01-16","USA/CUIMC-NP-5364/2021|MZ680974.1|2021-02-15","England/LOND-1359F9C/2021|OV314953.1|2021-11-25","USA/KPL-313_Nasal_Swab/2021|OP222737.1|2021-12-31","IMS-10150-CVDP-E66912C7-5ACA-4977-82E2-F6E1A70B5626|OU156786.1|2021-04-27","England/PLYM-20729CC/2021|OU864667.1|2021-10-06","node_1298822","Wales/PHWC-PDC7YQ/2021|2021-08-30","USA/NE-Noblis-S964.BC78.BT22.v1.1/2021|OM938984.1|2021-12-01","Denmark/DCGC-428071/2022|OX066483.1|2022-03-10","IMS-10183-CVDP-CADB58B9-51E0-43B4-9C19-7929C70E0374|OU567904.1|2021-07-24","AUS/VIC10358/2020|MW155076.1|2020-09-01","IMS-10087-CVDP-BA589A66-531B-4195-8258-22F9CC7C4497|OV679107.1|2021-12-14","node_933420","USA/CA-CDPH-2000010508/2020|ON216856.1|2020-12-11","Switzerland/SH-ETHZ-34783124/2021|OV326774.1|2021-11-14","USA/IL-CDC-LC0734800/2022|ON941006.1|2022-06-21","Northern_Ireland/NIRE-017cc7/2021|2021-12-17","England/ALDP-1778CF6/2021|OU358113.1|2021-06-18","USA/FL-CDC-ASC210249230/2021|MZ995678.1|2021-07-27","England/DHSC-CYDBKF5/2022|2022-01-04","USA/CO-CDPHE-2102916618/2022|OM860194.1|2022-02-10","USA/CA-CDC-ASC210113480/2021|OK296679.1|2021-07-11","THA/CONI-3328/2022|OP060789.1|2022-03-27","node_1179775","IMS-10023-CVDP-8806CB6B-FCB4-415E-BA08-B698B100AA90|OU109880.1|2021-03-21","BRA/SARS-CoV-2_mutants_generated_from_sequential_passages_in_Huh-7_cells/2020|OQ048259.1|2020-08-10","England/DHSC-CYYE8FU/2022|2022-02-20","England/PHEC-YYD1DP5/2022|2022-03-11","USA/CO-CDC-ASC210430719/2021|OL927628.1|2021-11-22","Switzerland/BL-ETHZ-35963862/2022|OV882561.1|2022-01-24","node_163491","England/QEUH-135D639/2021|OU034669.1|2021-02-28","node_880064","England/CAMC-C98EAD/2020|OD917787.1|2020-12-17","USA/WA-CDC-UW21121709950/2021|OM065614.1|2021-12-17","20220201815|OX448615.1|2022-01-11","OX314535.1|2022-01-01","USA/WI-UW-11513/2022|OP515931.1|2022-05-31","USA/SEARCH-0247-SAN/2020|MT810954.1|2020-03-27","USA/JW6938/2022|ON953779.1|2022-06-22","USA/FL-CDC-ASC210639430/2021|OM173419.1|2021-12-22","node_572285","hCoV-19/USA/NY-MSHSPSP-PV24650/2020|EPI_ISL_1300881|2020-12-06","USA/OK-CDC-ASC210638269/2021|OM172730.1|2021-12-20","USA/WI-CDC-ASC210638273/2021|OM172744.1|2021-12-21","node_691758","IMS-10172-CVDP-6B55DDB4-06C6-4887-A529-748794394385|OU122101.1|2021-03-10","node_357298","USA/FL-CDC-ASC210642292/2021|OM199050.1|2021-12-23","RNA|OW122540.1|2022-02-15","Molecular_surveillance_of_SARS-CoV-2_in_Germany|OV356720.1|2021-09-20","node_927614","USA/FL-BPHL-10749/2021|OL344995.1|2021-06-10","USA/NJ-CDC-LC0802670/2022|OP167632.1|2022-07-25","USA/IN-GD-083121-371/2021|OK158676.1|2021-08-14","hCoV-19/England/MILK-2DF642C/2021|EPI_ISL_7718520|2021-12-09","UZB/hCoV-19-Uzb-CGB-38/2021|MZ892624.1|2021-07-23","England/MILK-360C621/2022|OV969049.1|2022-02-08","node_1004504","USA/MI-CDC-STM-9UQJNQY2H/2021|OL499938.1|2021-10-21","node_911967","USA/ID-CDC-ASC210716187/2021|OL945444.1|2021-11-27","USA/SC-CDC-LC0013080/2021|MW642988.1|2021-02-01","USA/SC-CDC-ASC210244598/2021|MZ981235.1|2021-07-24","USA/MN-CDC-ASC210125244/2021|OK229736.1|2021-09-03","England/MILK-15268AA/2021|OU260439.1|2021-04-18","node_1250243","USA/NY-PRL-2021_0429_01C14/2021|MZ641224.1|2021-04-27","IRQ/Samawa-29/2021|MZ145315.1|2021-02-20","England/PHEP-020369/2021|2021-07-28","node_275296","USA/IL-CDC-ASC210584547/2022|OM517283.1|2022-01-13","Germany/2022|OW656430.1|2022-04-10","USA/CO-Noblis-S686B18/2021|MZ842997.1|2021-06-21","England/PHEC-3K044KB2/2021|2021-09-30","England/LSPA-3752852/2022|OW057126.1|2022-02-16","HUN/Hun-1/2020|OM812693.1|2020-11","node_905025","USA/CA-CDPH-FS48074591/2022|OQ112266.1|2022-11-16","USA/TG920857/2020|ON415701.1|2020-08-26","USA/TX-CDC-ASC210567824/2021|OM156012.1|2021-12-22","England/DHSC-CYBG3DE/2022|2022-02-04","USA/MN-CDC-IBX003537889913/2022|OM722843.1|2022-02-03","node_1127092","England/PHEC-4G0BAZ31/2021|2021-12-27","England/PHEC-4X03EZ00/2022|2022-01-22","USA/22007788/2022|OM614677.1|2022-01-19","England/DHSC-CYY4RZC/2022|2022-02-15","England/DHSC-CYYGBY8/2022|2022-01-04","England/MILK-32A2269/2022|OV734831.1|2022-01-15","USA/LA-EVTL9564/2021|OM997223.1|2021-12-13","USA/DE-CDC-ASC210560083/2021|OM331912.1|2021-12-28","USA/FL-CDC-ASC210571376/2021|OM261674.1|2021-12-27","OX307087.1|2021-11-26","USA/FL-BPHL-1553/2020|MW056104.1|2020-09-09","England/PHEC-YY84T8F/2023|2023-01-07","Molecular_surveillance_of_SARS-CoV-2_in_Germany|OV196096.1|2021-10-18","England/DHSC-CYN6EGB/2022|2022-01-16",
            "hCoV-19/England/MILK-9E05B3/2020|EPI_ISL_601443|2020-09-20","hCoV-19/Japan/IC-0564/2021|EPI_ISL_792683|2021-01-02","hCoV-19/USA/NY-MSHSPSP-PV24650/2020|EPI_ISL_1300881|2020-12-06","hCoV-19/India/MH-NCCS-P1162000182735/2021|EPI_ISL_1544014|2021-02-27","hCoV-19/Hong","hCoV-19/Australia/QLD2568/2021|EPI_ISL_7190366|2021-12-01","hCoV-19/England/MILK-2DF642C/2021|EPI_ISL_7718520|2021-12-09","hCoV-19/Denmark/DCGC-493190/2022|EPI_ISL_12248637|2022-04-17"
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
    }
    
    void analyze() {
        for (;!step();) {
            printf("\n");
        }

        // this->postprocess();

        // without neighbors
        this->em_postprocess(false);
        this->em_postprocess(true);
        // std::set<AuxNode*, MutationComparator> prev_iteration;
        // do {
        //     prev_iteration = this->selected_peaks;
        //     // with neighbors
        //     this->em_postprocess(true);
        // } while (prev_iteration != this->selected_peaks);
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
    /*
    peak_nodes.reserve(peak_nodes.size() + neighbor_nodes.size());
    peak_nodes.insert(peak_nodes.end(), neighbor_nodes.begin(), neighbor_nodes.end());
    neighbor_nodes.clear();
    
    generateFilteringData(T, T_condensed, condensed_node_mappings, ref_seq, peak_nodes, read_map, barcode_file, read_mutation_depth_vcf, condensed_nodes_csv);
    condensed_node_mappings.clear();
    MAT::clear_tree(T_condensed);
    */
}