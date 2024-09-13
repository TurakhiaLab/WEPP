#include "arena.hpp"

haplotype *arena::from_mat(haplotype *parent, MAT::Node *node)
{
    this->nodes.emplace_back();
    haplotype *ret = &this->nodes.back();
    ret->parent = parent;
    ret->mapped = false;
    ret->is_leaf = condensed_node_mappings.at(node).front()->is_leaf();
    ret->score = 0;
    ret->dist_divergence = 1;
    ret->stack_muts = {};
    ret->id = node->identifier;
    ret->condensed_source = node;
    ret->leaf_count = get_num_leaves(condensed_node_mappings, node);
    /* flatten mutation list */
    if (parent)
    {
        for (const MAT::Mutation &mut : parent->stack_muts)
        {
            /* only need to compare to our muts since parent's are unique */
            bool valid = true;
            for (const MAT::Mutation &comp : node->mutations)
            {
                if (comp.position == mut.position)
                {
                    valid = false;
                    break;
                }
            }

            if (valid)
                ret->stack_muts.push_back(mut);
        }
    }
    /* maybe rewinded mutation back to original */
    for (size_t i = 0; i < node->mutations.size(); ++i)
    {
        if (node->mutations[i].ref_nuc != node->mutations[i].mut_nuc)
        {
            ret->stack_muts.push_back(node->mutations[i]);
        }
    } 
    ret->muts = node->mutations;
    std::sort(ret->stack_muts.begin(), ret->stack_muts.end());
    std::sort(ret->muts.begin(), ret->muts.end());

    for (MAT::Node *child : node->children)
    {
        haplotype *curr = this->from_mat(ret, child);
        ret->children.emplace_back(curr);
    }

    return ret;
}
int arena::mat_tree_size(MAT::Node *node)
{
    int ret = 1;
    for (MAT::Node *child : node->children)
    {
        ret += mat_tree_size(child);
    }
    return ret;
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

void arena::print_cooccuring_mutations(int window_size)
{
    std::unordered_set<int> site_read_map;
    for (size_t i = 0; i < this->raw_reads.size(); i++)
    {
        const auto &rp = raw_reads[i];
        std::unordered_set<int> unknown_sites;
        for (const auto& mut: rp.mutations) {
            if (mut.mut_nuc == 0b1111)
                unknown_sites.insert(mut.position);
        }
        for (int j = rp.start; j <= rp.end; j++) {
            if (unknown_sites.find(j) == unknown_sites.end())
                site_read_map.insert(j);
        }
        unknown_sites.clear();
    }

    std::vector<std::string> comp = this->ds.true_haplotypes();
    for (const std::string &reference : comp)
    {
        auto sample_mutations = get_mutations(this->mat, reference);
        // Remove mutations from sample_mutations that are not present in site_read_map
        auto mut_itr = sample_mutations.begin();
        while (mut_itr != sample_mutations.end())
        {
            if (site_read_map.find(mut_itr->position) == site_read_map.end())
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
                if (sample_mutations[j].position - sample_mutations[i].position > window_size)
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
    std::unordered_set<int> site_read_map;
    for (size_t i = 0; i < this->raw_reads.size(); i++)
    {
        const auto &rp = raw_reads[i];
        std::unordered_set<int> unknown_sites;
        for (const auto& mut: rp.mutations) {
            if (mut.mut_nuc == 0b1111)
                unknown_sites.insert(mut.position);
        }
        for (int j = rp.start; j <= rp.end; j++) {
            if (unknown_sites.find(j) == unknown_sites.end())
                site_read_map.insert(j);
        }
        unknown_sites.clear();
    }

    double average_dist = 0.0;
    std::vector<std::string> comp = this->ds.true_haplotypes();
    printf("--- ground truth %zu haps; predicted %zu haps\n", comp.size(), selected.size());

    for (const std::string &reference : comp)
    {
        auto sample_mutations = get_mutations(this->mat, reference);
        // Remove mutations from sample_mutations that are not present in site_read_map
        auto mut_itr = sample_mutations.begin();
        while (mut_itr != sample_mutations.end())
        {
            if (site_read_map.find(mut_itr->position) == site_read_map.end())
                mut_itr = sample_mutations.erase(mut_itr);
            else
                mut_itr++;
        }

        MAT::Node *best_node = NULL;
        int min_dist = INT_MAX;
        for (const auto &pn : selected)
        {
            // Getting node_mutations from the Tree
            auto node_mutations = get_mutations(this->mat, condensed_node_mappings[pn->condensed_source].front()->identifier);
            // Remove mutations from node_mutations that are not present in site_read_map
            auto mut_itr = node_mutations.begin();
            while (mut_itr != node_mutations.end())
            {
                if (site_read_map.find(mut_itr->position) == site_read_map.end())
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
        printf("* dist: %02d true_node: %s pred_node: %s \n", min_dist, reference.c_str(), best_node->identifier.c_str());

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

void arena::print_flipped_mutation_distance(const std::vector<std::pair<haplotype *, double>> &selected)
{
    std::unordered_set<int> site_read_map;
    for (size_t i = 0; i < this->raw_reads.size(); i++)
    {
        const auto &rp = raw_reads[i];
        std::unordered_set<int> unknown_sites;
        for (const auto& mut: rp.mutations) {
            if (mut.mut_nuc == 0b1111)
                unknown_sites.insert(mut.position);
        }
        for (int j = rp.start; j <= rp.end; j++) {
            if (unknown_sites.find(j) == unknown_sites.end())
                site_read_map.insert(j);
        }
        unknown_sites.clear();
    }
    
    double average_dist = 0.0;
    std::vector<std::string> comp = this->ds.true_haplotypes();
    printf("--- ground truth %zu haps; predicted %zu haps\n", comp.size(), selected.size());

    for (const auto &pn : selected)
    {
        // Getting node_mutations from the Tree
        auto node_mutations = get_mutations(this->mat, pn.first->id);
        // Remove mutations from node_mutations that are not present in site_read_map
        auto mut_itr = node_mutations.begin();
        while (mut_itr != node_mutations.end())
        {
            if (site_read_map.find(mut_itr->position) == site_read_map.end())
                mut_itr = node_mutations.erase(mut_itr);
            else
                mut_itr++;
        }

        std::string best_node = "";
        int min_dist = INT_MAX;
        for (const std::string &reference : comp)
        {   
            std::vector<MAT::Mutation>sample_mutations;
            // Use reference tree for sample nodes if given as an input
            if (this->ref_mat.root != NULL)
                sample_mutations = get_mutations(this->ref_mat, reference);
            else
                sample_mutations = get_mutations(this->mat, reference);
            // Remove mutations from sample_mutations that are not present in site_read_map
            auto mut_itr = sample_mutations.begin();
            while (mut_itr != sample_mutations.end())
            {
                if (site_read_map.find(mut_itr->position) == site_read_map.end())
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

        std::string lineage_name = "";
        for (auto anc : mat.rsearch(condensed_node_mappings[pn.first->condensed_source].front()->identifier, true))
        {
            const auto &clade = anc->clade_annotations[1];
            if (clade != "")
            {
                lineage_name = clade;
                break;
            }
        }

        printf("PEAK: %s Lineage: %s weighted_dist: %.2f raw_dist: %02d proportion: %.2f (to) %s \n", pn.first->id.c_str(), lineage_name.c_str(), (min_dist * pn.second), min_dist, pn.second, best_node.c_str()); 
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
        std::string lineage_name;
        for (auto anc : mat.rsearch(condensed_node_mappings[p.first->condensed_source].front()->identifier, true))
        {
            const auto &clade = anc->clade_annotations[1];
            if (clade != "")
            {
                lineage_name = clade;
                break;
            }
        }
        a_map[lineage_name] += p.second;
    }

    std::cout << "--- lineage abundance " << std::endl;
    for (const auto &[name, val] : a_map)
    {
        printf("* lineage: %s abundance: %.6f\n", name.c_str(), val);
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