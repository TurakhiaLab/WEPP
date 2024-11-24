#include "arena.hpp"
#include "timer.hpp"

arena::arena(const dataset &ds) : ds{ds}, mat{ds.mat()}, coord{ds.mat()}
{
    // Masking mutations on reads
    this->masked_sites = ds.masked_sites();
    this->raw_reads = ds.reads();
    for (auto& rd: this->raw_reads)
    {
        auto mut_itr = rd.mutations.begin();
        while (mut_itr != rd.mutations.end())
        {
            if (std::find(this->masked_sites.begin(), masked_sites.end(), mut_itr->pos) != this->masked_sites.end())
                mut_itr = rd.mutations.erase(mut_itr);
            else
                mut_itr++;
        }
    }

    // (note that there is typically overallocation by condensation factor)
    // but shrink to fit may lead to pointer invalidation
    panmanUtils::Node* real_root = this->mat.root;
    this->nodes.reserve(pan_tree_size(real_root));
    std::vector<panmanUtils::Node*> empty;
    std::vector<mutation> empty_n_muts;

    timer t;
    this->from_pan(nullptr, real_root, this->site_read_map(), empty, empty_n_muts);
    std::cout << "Linearizing panmat took " << t.seconds() << " seconds " << std::endl << std::endl;
}

haplotype *
arena::from_pan(haplotype *parent, panmanUtils::Node *node, const std::unordered_set<int> &site_read_map, std::vector<panmanUtils::Node *> &parent_mapping, const std::vector<mutation>& condensed_n_muts)
{
    // if no mutations in site read map, condense and continue
    std::vector<mutation> muts;
    muts.reserve(node->nucMutation.size() + condensed_n_muts.size());
    this->get_single_mutations(muts, node);
    for (const auto& m: condensed_n_muts) 
    {
        auto it = std::lower_bound(muts.begin(), muts.end(), m);
        if (it != muts.end()) {
            if (it->pos != m.pos)
                muts.insert(it, m);
        }
        else {
            muts.emplace_back(m);
        }
    }
    tbb::parallel_sort(muts.begin(), muts.end());
    bool has_any = parent == nullptr; // root always gets added
    std::vector<mutation> child_condensed_n_muts;
    for (const mutation& mut: muts) {
        if (site_read_map.find(mut.pos) != site_read_map.end()) 
        {   
            if (mut.mut != NUC_N) {
                has_any = true;
                break;
            }
            else {
                child_condensed_n_muts.emplace_back(mut);
            }
        }
    }  
    
    if (!has_any) {
        parent_mapping.push_back(node);

        for (panmanUtils::Node *child : node->children)
        {
            haplotype *curr = this->from_pan(parent, child, site_read_map, parent_mapping, child_condensed_n_muts);
            if (curr) {
                parent->children.emplace_back(curr);
            }
        }
        return nullptr;
    }
    else
        child_condensed_n_muts.clear();

    this->nodes.emplace_back();
    haplotype *ret = &this->nodes.back();
    ret->depth = parent ? parent->depth + 1 : 0;
    ret->parent = parent;
    ret->mapped = false;
    ret->score = 0;
    ret->orig_score = 0;
    ret->dist_divergence = 1;
    ret->id = node->identifier;
    ret->condensed_source = node;
    ret->muts = {};
    ret->n_muts = {};
    ret->stack_muts = {};
    ret->stack_n_muts = {};
    ret->reference = &this->reference();

    // Update stack_muts
    if (parent)
    {
        // Update stack_muts from parent
        for (const mutation &mut : parent->stack_muts)
        {
            /* only need to compare to our muts since parent's are unique */
            bool valid = true;
            for (const mutation &comp : muts)
            {
                if (comp.pos == mut.pos)
                {
                    valid = false;
                    break;
                }
            }

            if (valid) {
                ret->stack_muts.push_back(mut);
            }
        }
        
        // Update stack_n_muts from parent
        for (const auto &coords : parent->stack_n_muts)
        {
            /* only need to compare to our muts since parent's are unique */
            mutation search;
            search.pos = coords.first;
            std::vector<int> break_mutations;
            auto it = std::lower_bound(muts.begin(), muts.end(), search);
            while ((it != muts.end()) && (it->pos <= coords.second)) 
            {
                break_mutations.emplace_back(it->pos);
                it++;
            }
            if (break_mutations.empty())
                ret->stack_n_muts.push_back(coords);
            else 
            {
                int start_point = coords.first, end_point;
                for (const auto& pos: break_mutations) 
                {
                    end_point = pos - 1;
                    ret->stack_n_muts.push_back({start_point, end_point});
                    start_point = pos + 1;
                }
                ret->stack_n_muts.push_back({start_point, coords.second});
            }
        }
    }
    
    /* maybe rewinded mutation back to original */
    for (const mutation &mut : muts)
    {
        if ((site_read_map.find(mut.pos) != site_read_map.end()) && (mut.ref != mut.mut))
        { 
            if (mut.mut != NUC_N)
            {
                ret->stack_muts.push_back(mut);
                ret->muts.push_back(mut);
            }
            else 
            {   
                // Update stack_n_muts
                auto it = std::lower_bound(ret->stack_n_muts.begin(), ret->stack_n_muts.end(), std::make_pair(mut.pos, INT_MIN),
                    [](const std::pair<int, int>& range, const std::pair<int, int>& value) {
                        return range.second < value.first; 
                });
                // No end coordinate is greater than mut.pos
                if (it == ret->stack_n_muts.end()) 
                {
                    if (it != ret->stack_n_muts.begin())
                    {
                        it = std::prev(it);
                        if ((mut.pos - it->second) == 1)
                            it->second = mut.pos;
                        else
                            ret->stack_n_muts.push_back({mut.pos, mut.pos}); 
                    }
                    else
                        ret->stack_n_muts.push_back({mut.pos, mut.pos}); 
                }
                else 
                {
                    // Check if can be merged at the start of selected block
                    if ((it->first - mut.pos) == 1)
                    {
                        it->first = mut.pos;
                        // Try merging with previous block
                        if (it != ret->stack_n_muts.begin()) {
                            auto prev_it = std::prev(it);
                            if ((mut.pos - prev_it->second) == 1) {
                                it->first = prev_it->first;
                                ret->stack_n_muts.erase(prev_it);
                            }
                        }
                    }
                    else 
                    {
                        // Check if can be merged at the end of previous block
                        if (it != ret->stack_n_muts.begin()) {
                            auto prev_it = std::prev(it);
                            if ((mut.pos - prev_it->second) == 1)
                                prev_it->second = mut.pos;
                            else
                                ret->stack_n_muts.insert(it, {mut.pos, mut.pos}); 
                        }
                        // Create a new block at the start of stack_n_muts
                        else
                            ret->stack_n_muts.insert(ret->stack_n_muts.begin(), {mut.pos, mut.pos});
                    }
                }

                // Update n_muts
                it = std::lower_bound(ret->n_muts.begin(), ret->n_muts.end(), std::make_pair(mut.pos, INT_MIN),
                    [](const std::pair<int, int>& range, const std::pair<int, int>& value) {
                        return range.second < value.first; 
                });
                // No end coordinate is greater than mut.pos
                if (it == ret->n_muts.end()) 
                {
                    if (it != ret->n_muts.begin())
                    {
                        it = std::prev(it);
                        if ((mut.pos - it->second) == 1)
                            it->second = mut.pos;
                        else
                            ret->n_muts.push_back({mut.pos, mut.pos}); 
                    }
                    else
                        ret->n_muts.push_back({mut.pos, mut.pos}); 
                }
                else 
                {
                    fprintf(stderr, "\n\nIncoming panmat mutations are NOT sorted\n\n");
                }
            }
        }
    }

    tbb::parallel_sort(ret->muts.begin(), ret->muts.end());
    tbb::parallel_sort(ret->n_muts.begin(), ret->n_muts.end());
    tbb::parallel_sort(ret->stack_muts.begin(), ret->stack_muts.end());
    tbb::parallel_sort(ret->stack_n_muts.begin(), ret->stack_n_muts.end());

    condensed_node_mappings[ret] = {node};
    std::vector<panmanUtils::Node*> &condensed_list = condensed_node_mappings[ret];
    for (panmanUtils::Node *child : node->children)
    {
        haplotype *curr = this->from_pan(ret, child, site_read_map, condensed_list, child_condensed_n_muts);
        if (curr) {
            ret->children.emplace_back(curr);
        }
    }

    return ret;
}

int arena::pan_tree_size(panmanUtils::Node *node)
{
    int ret = 1;
    for (panmanUtils::Node *child : node->children)
    {
        ret += pan_tree_size(child);
    }
    return ret;
}

std::vector<mutation> arena::get_mutations(const panmanUtils::Node* node, bool replace_GAP, bool replace_N) {
    std::vector<mutation> sample_mutations;
    std::unordered_map<int, mutation> parent_n_mutations;
    // Checking all ancestors of a node
    for (const panmanUtils::Node* curr = node; curr; curr = curr->parent) {
        std::vector<mutation> node_mutations;
        node_mutations.reserve(curr->nucMutation.size());
        get_single_mutations(node_mutations, curr);
        for (const auto& mut: node_mutations) {
            auto sm_itr = std::lower_bound(sample_mutations.begin(), sample_mutations.end(), mut);
            if (sm_itr == sample_mutations.end()) {
                sample_mutations.emplace_back(mut);
            }
            else if (sm_itr->pos != mut.pos) {
                sample_mutations.insert(sm_itr, mut);
            }
            else {
                if (((replace_GAP && (sm_itr->mut == NUC_GAP)) || (replace_N && (sm_itr->mut == NUC_N))) && (((!replace_GAP) || (mut.mut != NUC_GAP)) && ((!replace_N) || (mut.mut != NUC_N)))) {
                    if (parent_n_mutations.find(mut.pos) == parent_n_mutations.end())
                        parent_n_mutations[mut.pos] = mut;
                }
            }
        }
    }
    // Replace Gaps and remove back-mutations
    auto sm_itr = sample_mutations.begin();
    while (sm_itr != sample_mutations.end()) {
        //Update Gaps with their parent mutations
        if ((replace_GAP && (sm_itr->mut == NUC_GAP)) || (replace_N && (sm_itr->mut == NUC_N))) {
            if (parent_n_mutations.find(sm_itr->pos) != parent_n_mutations.end()) {
                const auto& par_mut = parent_n_mutations[sm_itr->pos];
                if (par_mut.ref != par_mut.mut) {
                    sm_itr->mut = par_mut.mut;
                    sm_itr++;
                }
                else {
                    sm_itr = sample_mutations.erase(sm_itr);
                }
            }
            else {
                sm_itr = sample_mutations.erase(sm_itr);
            }
        }
        // Remove Back-Mutations
        else if (sm_itr->ref == sm_itr->mut) {
            sm_itr = sample_mutations.erase(sm_itr);
        }
        else {
            sm_itr++;
        }
    }

    tbb::parallel_sort(sample_mutations.begin(), sample_mutations.end());
    
    return sample_mutations;
}

void arena::get_single_mutations(std::vector<mutation>& mutations, const panmanUtils::Node* node) {
    ::get_single_mutations(mutations, this->reference(), node, coord);
    return;
}

int arena::build_range_tree(int parent, haplotype *curr, int start, int end)
{
    int ret;
    if (parent == -1 || curr->has_mutations_in_range(start, end))
    {
        ret = ranged_nodes.size();
        ranged_nodes.emplace_back();
        ranged_nodes[ret].root = curr;
        ranged_nodes[ret].sources = {curr};
    }
    else
    {
        ret = parent;
        ranged_nodes[parent].sources.push_back(curr);
    }

    for (haplotype *child : curr->children)
    {
        int const sub = build_range_tree(ret, child, start, end);
        if (sub != ret)
        {
            ranged_nodes[ret].children.emplace_back(sub);
        }
    }

    return ret;
}

void arena::build_range_trees()
{
    if (!ranged_nodes.empty())
    {
        return;
    }

    // get ranges (heuristic based approach)
    using range = std::pair<int, int>;
    std::vector<range> read_ranges;
    std::transform(raw_reads.begin(), raw_reads.end(), std::back_inserter(read_ranges),
                   [](const raw_read &read)
                   {
                       return std::make_pair(read.start, read.end);
                   });
    std::sort(read_ranges.begin(), read_ranges.end());

    int const num = std::min((int)read_ranges.size(), NUM_RANGE_TREES);
    for (int i = 0; i < num; ++i)
    {
        int s = read_ranges.size() * i / num;
        int e = read_ranges.size() * (i + 1) / num;
        int read_start = read_ranges[s].first;
        int read_end = 0;
        for (int j = s; j < e; ++j)
        {
            read_end = std::max(read_end, read_ranges[j].second);
        }

        // there is technically a chance we can have a collision
        // where multiple range trees are for the exact same range
        // (to see this, if all reads are the exact same range
        // then all range trees will be the exact same)
        // However, this is not an issue even in the case that there are duplicates,
        // it only reduces performance, it doesn't actually make anything break
        range r{read_start, read_end};
        if (this->ranged_root_map.find(r) == ranged_root_map.end())
        {
            this->ranged_root_map[r] = build_range_tree(-1, &nodes[0], read_start, read_end);
        }
    }

    this->num_reads = 0;
    std::fill(true_read_counts.begin(), true_read_counts.end(), 0);

    read_distribution_bin_size = this->genome_size() / NUM_RANGE_BINS;
    for (const auto &r : raw_reads)
    {
        true_read_counts[r.start / read_distribution_bin_size] += r.degree;
        this->num_reads += r.degree;
    }
    for (size_t i = 0; i < NUM_RANGE_BINS; ++i)
    {
        true_read_distribution[i] = (double)true_read_counts[i] / this->num_reads;
    }
}

multi_haplotype *arena::find_range_tree_for(const raw_read &read)
{
    std::pair<int, int> search = {read.start, INT_MAX};
    auto it = ranged_root_map.upper_bound(search);
    do
    {
        it = std::prev(it);
        if (read.start >= it->first.first && read.end <= it->first.second)
        {
            return &ranged_nodes[it->second];
        }
    } while (it != ranged_root_map.begin());

    assert(0);
    return nullptr;
}

std::set<haplotype *, score_comparator> arena::closest_neighbors(haplotype *target, int max_radius, int num_limit) const
{
    std::set<haplotype *, score_comparator> all_neighbors, ret;
    std::queue<haplotype *> q;
    q.push(target);
    while (!q.empty())
    {
        haplotype *curr = q.front();
        q.pop();
        if (all_neighbors.find(curr) != all_neighbors.end() || curr->mutation_distance(target) > max_radius)
        {
            continue;
        }
        else
        {
            all_neighbors.insert(curr);
        }

        if (curr->parent)
        {
            q.push(curr->parent);
        }

        for (haplotype *hap : curr->children)
        {
            q.push(hap);
        }
    }

    // Find the num_limit-th element or the end of the set
    auto itr_end = (int)all_neighbors.size() > num_limit ? std::next(all_neighbors.begin(), num_limit) : all_neighbors.end();

    // Move the first num_limit elements to the ret set
    std::move(all_neighbors.begin(), itr_end, std::inserter(ret, ret.end()));

    return ret;
}

std::set<haplotype *, score_comparator> arena::highest_scoring_neighbors(haplotype *target, bool include_mapped, int max_radius, int num_limit) const
{
    auto dfs_possible_neighbors = [&](auto &&dfs, haplotype *pivot, haplotype *curr, std::set<haplotype *, score_comparator> &s)
    {
        if (pivot->mutation_distance(curr) > max_radius)
        {
            return;
        }

        if (!curr->mapped || include_mapped)
        {
            s.insert(curr);
        }

        for (haplotype *child : curr->children)
        {
            dfs(dfs, pivot, child, s);
        }
    };

    /* find all possible neighbors within a radius from a given node */
    auto possible_neighbors = [&](haplotype *pivot, std::set<haplotype *, score_comparator> &s)
    {
        haplotype *curr = pivot;
        while (curr->parent && pivot->mutation_distance(curr->parent) <= max_radius)
        {
            curr = curr->parent;
        }

        dfs_possible_neighbors(dfs_possible_neighbors, pivot, curr, s);
    };

    std::set<haplotype *, score_comparator> ret;
    possible_neighbors(target, ret);
    while ((int)ret.size() > num_limit)
    {
        ret.erase(std::prev(ret.end()));
    }

    return ret;
}

void arena::print_mutation_distance(const std::vector<haplotype *> &selected)
{
    std::unordered_set<int> site_read_map = this->site_read_map();

    double average_dist = 0.0;
    std::vector<std::string> comp = this->ds.true_haplotypes();
    printf("--- ground truth %zu haps; predicted %zu haps\n", comp.size(), selected.size());

    std::vector<std::vector<mutation>> selected_mutations;
    for (const auto &pn : selected) {
        auto node_mutations = get_mutations(pn->condensed_source);
        
        // Remove mutations from node_mutations that are not present in site_read_map
        auto mut_itr = node_mutations.begin();
        while (mut_itr != node_mutations.end())
        {
            if (site_read_map.find(mut_itr->pos) == site_read_map.end())
                mut_itr = node_mutations.erase(mut_itr);
            else
                mut_itr++;
        }
       selected_mutations.emplace_back(node_mutations);
    }

    for (const std::string &reference : comp)
    {
        std::vector<mutation>sample_mutations = get_mutations(this->mat.allNodes.at(reference));

        // Remove mutations from sample_mutations that are not present in site_read_map
        auto mut_itr = sample_mutations.begin();
        while (mut_itr != sample_mutations.end())
        {
            if (site_read_map.find(mut_itr->pos) == site_read_map.end())
                mut_itr = sample_mutations.erase(mut_itr);
            else
                mut_itr++;
        }

        std::string best_node = "";
        int min_dist = INT_MAX;
        for (size_t i = 0; i < selected.size(); i++)
        {
            const auto& pn = selected[i];
            auto node_mutations = selected_mutations[i];
            int curr_dist = mutation_distance(sample_mutations, node_mutations);
            if (curr_dist < min_dist)
            {
                min_dist = curr_dist;
                best_node = pn->id;
            }
        }

        average_dist += (double)min_dist / comp.size();
        printf("* dist: %02d true_node: %s pred_node: %s \n", min_dist, reference.c_str(), best_node.c_str());

    }

    printf("average mutation_distance: %0.3f\n", average_dist);
}

void arena::print_flipped_mutation_distance(const std::vector<std::pair<haplotype *, double>> &selected)
{
    std::unordered_set<int> site_read_map = this->site_read_map();

    double average_dist = 0.0;
    std::vector<std::string> comp = this->ds.true_haplotypes();
    printf("--- ground truth %zu haps; predicted %zu haps\n", comp.size(), selected.size());

    std::vector<std::vector<mutation>> comp_mutations;
    for (const std::string &reference : comp) {
        auto sample_mutations = get_mutations(this->mat.allNodes.at(reference));
        
        // Remove mutations from sample_mutations that are not present in site_read_map
        auto mut_itr = sample_mutations.begin();
        while (mut_itr != sample_mutations.end())
        {
            if (site_read_map.find(mut_itr->pos) == site_read_map.end())
                mut_itr = sample_mutations.erase(mut_itr);
            else
                mut_itr++;
        }
       comp_mutations.emplace_back(sample_mutations);
    }
    
    for (const auto &pn : selected)
    {
        // Getting node_mutations from the Tree
        auto node_mutations = get_mutations(pn.first->condensed_source);
        
        // Remove mutations from node_mutations that are not present in site_read_map
        auto mut_itr = node_mutations.begin();
        while (mut_itr != node_mutations.end())
        {
            if (site_read_map.find(mut_itr->pos) == site_read_map.end())
                mut_itr = node_mutations.erase(mut_itr);
            else
                mut_itr++;
        }

        std::string best_node = "";
        int min_dist = INT_MAX;
        for (size_t i = 0; i < comp.size(); i++)
        {
            auto reference = comp[i];
            auto sample_mutations = comp_mutations[i];
            int curr_dist = mutation_distance(sample_mutations, node_mutations);
            if (curr_dist < min_dist)
            {
                min_dist = curr_dist;
                best_node = reference;
            }
        }

        std::string lineage_name = "";
        for (auto anc = condensed_node_mappings[pn.first].front(); anc; anc = anc->parent)
        {
            if (anc->annotations.size()) {
                lineage_name = anc->annotations.front();
                break;
            }
        }

        printf("PEAK: %s Lineage: %.*s weighted_dist: %.2f raw_dist: %02d proportion: %.2f (to) %s \n", pn.first->id.c_str(), static_cast<int>(lineage_name.size() - 1), lineage_name.c_str(), (min_dist * pn.second), min_dist, pn.second, best_node.c_str()); 
        average_dist += (min_dist * pn.second);
    }

    printf("average (flipped) mutation_distance: %0.3f\n", average_dist);
}

void arena::print_full_report(const std::vector<std::pair<haplotype *, double>> &abundance)
{
    std::cout << "----- [final report] -----" << std::endl
              << std::endl;

    print_flipped_mutation_distance(abundance);

    std::unordered_map<std::string, double> a_map;
    for (const auto &p : abundance)
    {
        std::string lineage_name = "";
        for (auto anc = condensed_node_mappings[p.first].front(); anc; anc = anc->parent)
        {
            if (anc->annotations.size()) {
                lineage_name = anc->annotations.front();
                break;
            }
        }
        
        if (!lineage_name.empty()) {
            a_map[lineage_name] += p.second;
        }
    }

    std::cout << "--- lineage abundance " << std::endl;
    for (const auto &[name, val] : a_map)
    {
        printf("* lineage: %.*s abundance: %.6f\n", static_cast<int>(name.size() - 1), name.c_str(), val);
    }
}

void arena::dump_haplotype_proportion(const std::vector<std::pair<haplotype *, double>> &abundance)
{
    std::ofstream csv(this->ds.haplotype_proportion_path());
    std::string csv_print;
    for (const auto &n_p : abundance)
    {
        csv_print = n_p.first->id + "," + std::to_string(n_p.second);
        csv_print += "\n";
        csv << csv_print;
    }
}

void arena::dump_read2haplotype_mapping(const std::vector<std::pair<haplotype *, double>> &abundance)
{
    std::ofstream csv(this->ds.haplotype_read_path());
    std::string csv_print;

    const std::vector<raw_read> &reads = this->reads();
    std::unordered_map<std::string, std::vector<std::string>> reverse_merge = this->ds.read_reverse_merge();

    tbb::concurrent_hash_map<std::string, std::vector<std::string>> haplotype_reads_map;
    static tbb::affinity_partitioner ap;
    tbb::parallel_for(tbb::blocked_range<size_t>(0, reads.size()), [&](tbb::blocked_range<size_t> range)
    {
        for (size_t i = range.begin(); i < range.end(); ++i) {
            //Find EPPs for every read
            std::vector<haplotype*> epps;
            float min_dist = std::numeric_limits<float>::max();
            for (const auto& curr_node: abundance) {
                float curr_dist = curr_node.first->mutation_distance(reads[i]);
                if (curr_dist <= min_dist) {
                    if (curr_dist < min_dist) {
                        min_dist = curr_dist;
                        epps.clear();
                    }
                    epps.emplace_back(curr_node.first);
                }
            }
            
            //Account for all the reads while assigning to EPPs
            const auto& read_names = reverse_merge[reads[i].read];
            for (const auto& curr_node: epps) {
                tbb::concurrent_hash_map<std::string, std::vector<std::string>>::accessor ac;
                auto created = haplotype_reads_map.insert(ac, std::make_pair(curr_node->id, read_names));
                if (!created)
                    ac->second.insert(ac->second.end(), read_names.begin(), read_names.end());
                ac.release();
            }
        } 
    }, ap);

    // Write CSV
    for (const auto &n_r : haplotype_reads_map)
    {
        csv_print = n_r.first;
        for (const auto &r_name : n_r.second)
            csv_print += "," + r_name;
        csv_print += "\n";
        csv << csv_print;
    }

    // Run Python script to generate sam files
    std::string command = "python src/sam_generation.py " + ds.directory() + " " + ds.file_prefix();
    int result = std::system(command.c_str());
    if (result)
        fprintf(stderr, "\nCannot run sam_generation.py\n");
}

void arena::resolve_unaccounted_mutations(const std::vector<std::pair<haplotype *, double>> &abundance)
{
    std::ofstream csv_haplotypes(this->ds.mutation_haplotypes_path()), csv_reads(this->ds.mutation_reads_path());
    std::string csv_print_haplotypes, csv_print_reads;

    // Read residual_mutations
    std::string file_path = this->ds.directory() + "/residual_mutations.txt";
    std::ifstream file(file_path);  
    std::vector<std::tuple<int, char, float>> mutations; 

    if (file.is_open()) {
        std::string line;
        // Read each line and store it in the mutations
        while (std::getline(file, line)) {
            std::size_t comma_pos = line.find(',');
            std::string mutation_part = line.substr(0, comma_pos);
            float value = std::stof(line.substr(comma_pos + 1));
            int pos = std::stoi(mutation_part.substr(0, mutation_part.size() - 1));
            char nuc = mutation_part.back();
            mutations.emplace_back(std::make_tuple(pos, nuc, value));
        }
        fprintf(stderr, "Residual mutations: %ld\n", mutations.size());
        file.close();  
    } 
    else {
        std::cerr << "Unable to open file: " << file_path << std::endl;
    }

    // Mask positions on the reads covered by mutations
    tbb::concurrent_hash_map<std::string, std::vector<raw_read>> mutations_read_map;
    tbb::concurrent_hash_map<std::string, std::vector<haplotype *>> mutations_haplotype_map;
    const std::vector<raw_read> &reads = this->reads();
    std::unordered_map<std::string, std::vector<std::string>> reverse_merge = this->ds.read_reverse_merge();

    static tbb::affinity_partitioner ap;
    tbb::parallel_for(tbb::blocked_range<size_t>(0, reads.size()), [&](const tbb::blocked_range<size_t>& r) 
    {
        for (size_t i = r.begin(); i < r.end(); ++i) {
            raw_read rp = reads[i];
            // Get residual sites covered by rp
            std::vector<std::tuple<int, char, float>> residual_sites_covered;
            std::copy_if(mutations.begin(), mutations.end(), std::back_inserter(residual_sites_covered),
                     [rp](const std::tuple<int, char, float>& mutation) {
                         return std::get<0>(mutation) >= rp.start && std::get<0>(mutation) <= rp.end;
                     });

            // Get residual_mutations_covered covered by reads and convert them to 'N'
            std::vector<std::string> residual_mutations_covered;
            for (const auto& mut_tuple: residual_sites_covered) {
                std::string curr_residual_mut = std::to_string(std::get<0>(mut_tuple)) + std::get<1>(mut_tuple) + ":" + std::to_string(std::get<2>(mut_tuple));
                bool site_found = false;
                for (auto &mut: rp.mutations) {
                    if (mut.pos == std::get<0>(mut_tuple)) {
                        // If NOT 'N' then both position and allele should match
                        if (mut.mut != NUC_N) {
                            if (char_from_nuc(mut.mut) == std::get<1>(mut_tuple)) {
                                // Masking the mutation on site before pushing
                                mut.mut = NUC_N;
                                residual_mutations_covered.emplace_back(curr_residual_mut);
                            }
                        }
                        // If 'N' then only position should match
                        else
                            residual_mutations_covered.emplace_back(curr_residual_mut);
                        site_found = true;
                        break;
                    }
                }
                if (!site_found) {
                    if (ds.reference()[std::get<0>(mut_tuple) - 1] == std::get<1>(mut_tuple)) {
                        // Add site as masked mutation
                        mutation mut;
                        mut.ref = nuc_from_char(std::get<1>(mut_tuple));
                        mut.mut = NUC_N;
                        mut.pos = std::get<0>(mut_tuple);
                        rp.mutations.emplace_back(mut);
                        std::sort(rp.mutations.begin(), rp.mutations.end());
                        residual_mutations_covered.emplace_back(curr_residual_mut);
                    }
                }
            }

            // Add rp to mutations_read_map if it covers residual mutations
            for (const auto& curr_residual_mut: residual_mutations_covered) {
                tbb::concurrent_hash_map<std::string, std::vector<raw_read>>::accessor ac;
                mutations_read_map.insert(ac, curr_residual_mut);
                ac->second.emplace_back(rp);
                ac.release();
            }
        }
    }, ap);
    
    // Write csv_reads
    for (const auto &m_r : mutations_read_map)
    {
        csv_print_reads = m_r.first;
        for (const auto &rp : m_r.second) {
            for (const auto& r_name: reverse_merge[rp.read])
                csv_print_reads += "," + r_name;
        }
        csv_print_reads += "\n";
        csv_reads << csv_print_reads;
    }

    // Find EPPs for each residual mutation in mutations_read_map
    tbb::parallel_for(tbb::blocked_range<size_t>(0, mutations_read_map.size()), [&](const tbb::blocked_range<size_t>& r) 
    {
        for (size_t i = r.begin(); i != r.end(); ++i) {
            tbb::concurrent_hash_map<std::string, std::vector<raw_read>>::const_accessor k_ac;
            auto itr = mutations_read_map.begin();
            std::advance(itr, i);

            if (mutations_read_map.find(k_ac, itr->first)) {  
                std::unordered_map<haplotype *, int> hap_reads_count;
                for (const auto& rp: k_ac->second) 
                {
                    //Find EPPs for every read
                    std::vector<haplotype* > epps;
                    float min_dist = std::numeric_limits<float>::max();
                    for (const auto& curr_node_abun: abundance) {
                        float curr_dist = curr_node_abun.first->mutation_distance(rp);
                        if (curr_dist <= min_dist) {
                            if (curr_dist < min_dist) {
                                min_dist = curr_dist;
                                epps.clear();
                            }
                            epps.emplace_back(curr_node_abun.first);
                        }
                    }

                    // Add read count to each EPP in hap_reads_count
                    for (const auto& curr_node: epps) {
                        if (hap_reads_count.find(curr_node) == hap_reads_count.end())
                            hap_reads_count[curr_node] = rp.degree;
                        else
                            hap_reads_count[curr_node] += rp.degree;
                    }
                }

                // Find max reads mapping to a haplotype
                int max_reads = 0;
                auto max_itr = std::max_element(hap_reads_count.begin(), hap_reads_count.end(),
                    [](const auto& lhs, const auto& rhs) {
                        return lhs.second < rhs.second;
                    });
                if (max_itr != hap_reads_count.end())
                    max_reads = max_itr->second;
                else
                    fprintf(stderr, "\nThere are ZERO reads mapping to the selected peaks\n\n");

                // Store unseen mutations and their haplotypes in mutations_haplotype_map
                std::vector<haplotype *> possible_haplotypes;
                for (const auto& h_c: hap_reads_count) {
                    if (h_c.second == max_reads)
                        possible_haplotypes.emplace_back(h_c.first);
                }

                // Add haplotypes for current mutation to mutations_haplotype_map
                tbb::concurrent_hash_map<std::string, std::vector<haplotype *>>::accessor ac;
                mutations_haplotype_map.insert(ac, std::make_pair(k_ac->first, possible_haplotypes));
                ac.release();
            }
        }
    }, ap);
    
    // Write csv_haplotypes
    for (const auto &m_h : mutations_haplotype_map)
    {
        csv_print_haplotypes = m_h.first;
        for (const auto &hap : m_h.second) {
            csv_print_haplotypes += "," + hap->id;
        }
        csv_print_haplotypes += "\n";
        csv_haplotypes << csv_print_haplotypes;
    }
}