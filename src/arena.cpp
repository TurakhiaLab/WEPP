#include "arena.hpp"
#include "timer.hpp"

arena::arena(const dataset &ds) : ds{ds}, mat{ds.mat()}, coord{ds.mat()}
{
    this->raw_reads = ds.reads();
    
    // (note that there is typically overallocation by condensation factor)
    // but shrink to fit may lead to pointer invalidation
    panmanUtils::Node* real_root = this->mat.root;
    this->nodes.reserve(pan_tree_size(real_root));
    std::vector<panmanUtils::Node*> empty;

    timer t;
    this->from_pan(nullptr, real_root, this->site_read_map(), empty);
    std::cout << "From pan took " << t.seconds() << " seconds " << std::endl;
}


haplotype *
arena::from_pan(haplotype *parent, panmanUtils::Node *node, const std::unordered_set<int> &site_read_map, std::vector<panmanUtils::Node *> &parent_mapping)
{
    // if no mutations in site read map, condense and continue
    std::vector<mutation> muts = this->get_single_mutations(node);
    bool has_any = parent == nullptr; // root always gets added
    for (const mutation& mut: muts) {
        if (site_read_map.find(mut.pos) != site_read_map.end()) {
            has_any = true;
            break;
        }
    }  

    if (!has_any) {
        parent_mapping.push_back(node);

        for (panmanUtils::Node *child : node->children)
        {
            haplotype *curr = this->from_pan(parent, child, site_read_map, parent_mapping);
            if (curr) {
                parent->children.emplace_back(curr);
            }
        }
        return nullptr;
    }

    this->nodes.emplace_back();
    haplotype *ret = &this->nodes.back();
    ret->depth = parent ? parent->depth + 1 : 0;
    ret->parent = parent;
    ret->mapped = false;
    ret->score = 0;
    ret->dist_divergence = 1;
    ret->id = node->identifier;
    ret->condensed_source = node;

    ret->muts = parent ? std::move(muts) : std::vector<mutation>{};
    std::sort(ret->muts.begin(), ret->muts.end());

    condensed_node_mappings[node] = {node};
    std::vector<panmanUtils::Node*> &condensed_map = condensed_node_mappings[node];
    for (panmanUtils::Node *child : node->children)
    {
        haplotype *curr = this->from_pan(ret, child, site_read_map, condensed_map);
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

std::vector<mutation> arena::get_mutations(const panmanUtils::Node* node) {
    std::vector<mutation> sample_mutations;
    for (const panmanUtils::Node* curr = node; curr; curr = curr->parent) { //Checking all ancestors of a node
        for (const auto& mut: get_single_mutations(curr)) {
            auto sm_itr = sample_mutations.begin();
            while (sm_itr != sample_mutations.end()) {
                if (sm_itr->pos == mut.pos)
                    break;
                sm_itr++;
            }
            if (sm_itr == sample_mutations.end()) {
                sample_mutations.emplace_back(mut);
            }
        }
    }
    //Remove Back-Mutations
    auto sm_itr = sample_mutations.begin();
    while (sm_itr != sample_mutations.end()) {
        if (sm_itr->ref == sm_itr->mut)
            sm_itr = sample_mutations.erase(sm_itr);
        else {
            sm_itr++;
        }
    }

    sort(sample_mutations.begin(), sample_mutations.end());
    
    return sample_mutations;
}

std::vector<mutation> arena::get_single_mutations(const panmanUtils::Node* node) {
    return ::get_single_mutations(this->reference(), node, coord);
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

std::set<haplotype *> arena::closest_neighbors(haplotype *target, int max_radius, int num_limit) const
{
    std::set<haplotype *> ret;

    std::queue<haplotype *> q;
    q.push(target);
    while (!q.empty())
    {
        haplotype *curr = q.front();
        q.pop();
        if (ret.find(curr) != ret.end() || curr->mutation_distance(target) > max_radius)
        {
            continue;
        }
        else
        {
            ret.insert(curr);
            if ((int)ret.size() >= num_limit)
            {
                return ret;
            }
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

void arena::print_cooccuring_mutations(int window_size)
{
    std::unordered_set<int> site_read_map = this->site_read_map();

    std::vector<std::string> comp = this->ds.true_haplotypes();
    for (const std::string &reference : comp)
    {
        auto sample_mutations = this->get_mutations(this->mat.allNodes.at(reference));
        // Remove mutations from sample_mutations that are not present in site_read_map
        auto mut_itr = sample_mutations.begin();
        while (mut_itr != sample_mutations.end())
        {
            if (site_read_map.find(mut_itr->pos) == site_read_map.end())
                mut_itr = sample_mutations.erase(mut_itr);
            else
                mut_itr++;
        }

        std::sort(sample_mutations.begin(), sample_mutations.end());

        std::cout << " reference: " << reference << std::endl;
        for (size_t i = 0; i < sample_mutations.size(); ++i)
        {
            for (size_t j = i + 1; j < sample_mutations.size(); ++j)
            {
                if (sample_mutations[j].pos - sample_mutations[i].pos > window_size)
                {
                    break;
                }

                std::cout << " first " << sample_mutations[i].get_string()
                          << " second << " << sample_mutations[j].get_string() << std::endl;
            }
        }
    }
}

void arena::print_mutation_distance(const std::vector<haplotype *> &selected)
{
    std::unordered_set<int> site_read_map = this->site_read_map();

    double average_dist = 0.0;
    std::vector<std::string> comp = this->ds.true_haplotypes();
    printf("--- ground truth %zu haps; predicted %zu haps\n", comp.size(), selected.size());

    for (const std::string &reference : comp)
    {
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

        panmanUtils::Node *best_node = NULL;
        int min_dist = INT_MAX;
        for (const auto &pn : selected)
        {
            auto node_mutations = get_mutations(condensed_node_mappings[pn->condensed_source].front());
            auto mut_itr = node_mutations.begin();
            while (mut_itr != node_mutations.end())
            {
                if (site_read_map.find(mut_itr->pos) == site_read_map.end())
                    mut_itr = node_mutations.erase(mut_itr);
                else
                    mut_itr++;
            }

            int curr_dist = mutation_distance(sample_mutations, node_mutations);
            if (curr_dist < min_dist)
            {
                min_dist = curr_dist;
                best_node = pn->condensed_source;
            }
        }

        average_dist += (double)min_dist / comp.size();
        printf("* dist: %02d true node: %s pred node: %s \n", min_dist, reference.c_str(), best_node->identifier.c_str());

        // auto node_mutations = get_mutations(this->mat, best_node->identifier);
        // // Remove mutations from node_mutations that are not present in site_read_map
        // mut_itr = node_mutations.begin();
        // while (mut_itr != node_mutations.end())
        // {
        //     if (site_read_map.find(mut_itr->position) == site_read_map.end())
        //         mut_itr = node_mutations.erase(mut_itr);
        //     else
        //         mut_itr++;
        // }
        // std::vector<MAT::Mutation> diffs;
        // std::sort(sample_mutations.begin(), sample_mutations.end());
        // std::sort(node_mutations.begin(), node_mutations.end());
        // std::set_symmetric_difference(sample_mutations.begin(), sample_mutations.end(), node_mutations.begin(), node_mutations.end(),
        //     std::back_inserter(diffs)
        // );
        // for (size_t i = 0; i < diffs.size(); ++i) {
        //     std::cerr << " Diff at " << diffs[i].get_string() << std::endl;
        // }
    }

    printf("average mutation_distance: %0.3f\n", average_dist);
}

void arena::print_flipped_mutation_distance(const std::vector<haplotype *> &selected)
{
    std::unordered_set<int> site_read_map = this->site_read_map();

    double average_dist = 0.0;
    std::vector<std::string> comp = this->ds.true_haplotypes();
    printf("--- ground truth %zu haps; predicted %zu haps\n", comp.size(), selected.size());

    for (const auto &pn : selected)
    {
        // Getting node_mutations from the Tree
        auto node_mutations = get_mutations(condensed_node_mappings[pn->condensed_source].front());
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
        for (const std::string &reference : comp)
        {
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

            int curr_dist = mutation_distance(sample_mutations, node_mutations);
            if (curr_dist < min_dist)
            {
                min_dist = curr_dist;
                best_node = reference;
            }
        }

        average_dist += (double)min_dist / selected.size();
        printf("* dist: %02d (to) %s \n", min_dist, best_node.c_str()); 
    }

    printf("average (flipped) mutation_distance: %0.3f\n", average_dist);
}


void arena::print_full_report(const std::vector<std::pair<haplotype *, double>> &abundance)
{
    std::cout << "----- [final report] -----" << std::endl
              << std::endl;

    std::vector<haplotype *> haps;
    std::transform(abundance.begin(), abundance.end(), std::back_inserter(haps),
                   [](const auto &p)
                   {
                       return p.first;
                   });
    
    print_mutation_distance(haps);

    std::unordered_map<std::string, double> a_map;
    for (const auto &p : abundance)
    {
        std::string lineage_name;
        panmanUtils::Node* target = condensed_node_mappings[p.first->condensed_source].front();
        for (panmanUtils::Node* curr = target; target; target = target->parent) {
            const auto &clade = curr->annotations[1];
            if (clade != "")
            {
                lineage_name = clade;
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
        printf("* lineage: %s abundance: %.6f\n", name.c_str(), val);
    }
}

void arena::dump_read2node_mapping(const std::vector<std::pair<haplotype *, double>> &abundance)
{
    std::ofstream csv(this->ds.haplotype_read_path());
    std::string csv_print;

    const std::vector<raw_read> &reads = this->reads();
    std::unordered_map<std::string, std::vector<std::string>> reverse_merge = ds.read_reverse_merge();

    tbb::concurrent_hash_map<std::string, std::vector<std::string>> node_reads_map;
    static tbb::affinity_partitioner ap;
    tbb::parallel_for(tbb::blocked_range<size_t>(0, reads.size()), [&](tbb::blocked_range<size_t> range)
                      {
            for (size_t i = range.begin(); i < range.end(); ++i) {
                //Find EPPs for every read
                std::vector<haplotype*> epps;
                int min_dist = std::numeric_limits<int>::max();
                for (const auto& curr_node: abundance) {
                    int curr_dist = curr_node.first->mutation_distance(reads[i]);
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
                    auto created = node_reads_map.insert(ac, std::make_pair(curr_node->id, read_names));
                    if (!created)
                        ac->second.insert(ac->second.end(), read_names.begin(), read_names.end());
                    ac.release();
                }
            } }, ap);

    for (const auto &n_r : node_reads_map)
    {
        csv_print = n_r.first;
        for (const auto &r_name : n_r.second)
            csv_print += "," + r_name;
        csv_print += "\n";
        csv << csv_print;
    }

    // Run Python script to generate sam files
    std::string command = "python src/WBE/sam_generation.py " + ds.directory() + " " + ds.file_prefix();
    int result = std::system(command.c_str());
    if (result)
        fprintf(stderr, "\nCannot run sam_generation.py\n");
}