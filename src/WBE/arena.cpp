#include "arena.hpp"

haplotype *arena::from_mat(haplotype *parent, MAT::Node *node)
{
    this->nodes.emplace_back();
    haplotype *ret = &this->nodes.back();
    ret->parent = parent;
    ret->mapped = false;
    ret->is_leaf = condensed_node_mappings.at(node).front()->is_leaf();
    ret->orig_score = 0;
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

void arena::get_residual_cooccuring_mutations(int window_size)
{
    // Read residual_mutations
    std::string file_path = this->ds.directory() + "/residual_mutations.txt";
    std::ifstream file(file_path);  
    std::vector<std::string> mutations;  

    if (file.is_open()) {
        std::string line;
        // Read each line and store it in the vector
        while (std::getline(file, line)) {
            mutations.emplace_back(line);
        }
        file.close();  
    } 
    else
        std::cerr << "Unable to open file: " << file_path << std::endl;

    // Sort mutations based on site
    std::sort(mutations.begin(), mutations.end(), [](const std::string &a, const std::string &b) {
        int a_pos = std::stoi(a.substr(0, a.size() - 1));
        int b_pos = std::stoi(b.substr(0, b.size() - 1));

        if (a_pos != b_pos)
            return std::stoi(a) < std::stoi(b);
        else 
            return a.back() < b.back();
    });

    // Create sets of co-occuring mutations
    std::vector<std::vector<std::string>> cooccuring_mutations;
    for (size_t i = 0; i < mutations.size(); i++)
    {
        std::vector<std::string> cooccuring;
        auto ref_mut = mutations[i];
        int ref_position = std::stoi(ref_mut.substr(0, ref_mut.size() - 1));
        for (size_t j = i + 1; j < mutations.size(); j++) 
        {
            auto cmp_mut = mutations[j];
            int cmp_position = std::stoi(cmp_mut.substr(0, cmp_mut.size() - 1));
            if ((cmp_position - ref_position) > window_size)
                break;
            else
                cooccuring.emplace_back(cmp_mut);
        }
        if (!cooccuring.empty()) 
        {
            cooccuring.emplace(cooccuring.begin(), ref_mut);
            cooccuring_mutations.emplace_back(cooccuring);
        }
    }

    // Check freq of co-occuring mutations
    for (auto cooccuring: cooccuring_mutations)
    {
        int start_position = std::stoi(cooccuring.front().substr(0, cooccuring.front().size() - 1));
        int end_position = std::stoi(cooccuring.back().substr(0, cooccuring.back().size() - 1));

        std::string grouped_mutations;
        std::vector<int> ref_sites;
        std::vector<std::string> mutations;
        for (auto mut: cooccuring) {
            if (mut != cooccuring.front())
                grouped_mutations += ":";
            grouped_mutations += mut;
            int site = std::stoi(mut.substr(0, mut.size() - 1));
            if (ds.reference()[site - 1] == mut.back())
                ref_sites.emplace_back(site);
            else
                mutations.emplace_back(mut);
        }

        // Checking reads
        int reads_mutations_present = 0, reads_sites_present = 0;
        std::vector<raw_read> reads_carrying_muts;
        for (size_t i = 0; i < this->raw_reads.size(); i++)
        {
            const auto &rp = raw_reads[i];
            if ((rp.start <= start_position) && (rp.end >= end_position)) {
                reads_sites_present += rp.degree;
                int mutations_covered = 0;
                for (auto mut: rp.mutations) {
                    if (mut.mut_nuc != 0b1111) {
                        auto ref_sites_itr = std::find(ref_sites.begin(), ref_sites.end(), mut.position);
                        if (ref_sites_itr != ref_sites.end())
                        {
                            mutations_covered = 0;
                            break;
                        }
                        else 
                        {
                            std::string curr_mut = std::to_string(mut.position) + MAT::get_nuc(mut.mut_nuc);
                            auto mutations_itr = std::find(mutations.begin(), mutations.end(), curr_mut);
                            if (mutations_itr != mutations.end())
                                mutations_covered++;
                        }
                    }
                }
                if (mutations_covered == (int)mutations.size()) {
                    reads_mutations_present += rp.degree;
                    reads_carrying_muts.emplace_back(rp);
                }
            }
        }

        if (reads_sites_present) 
        {
            printf("%s -> %f,%d,%d\n", grouped_mutations.c_str(), (double)reads_mutations_present / reads_sites_present, reads_mutations_present, reads_sites_present);

            std::vector<std::string> selected_nodes = {
			    	"node_991589",
	                "USA/NY-CDC-LC0978715/2022|OQ206817.1|2022-12-31",
	                "CHN/SH-XG2307-9897/2023|C_AA029986.1|2023-06-26",
	                "node_987817",
	                "England/PHEC-YYE5AXM/2023|2023-06-17",
	                "USA/IL-RIPHL_121818_G/2023|OR103062.1|2023-05-24",
	                "node_994586"
            };

            std::vector<haplotype> selected;
            std::vector<std::string> single_epps;
            for (auto hap: this->haplotypes())
            {
                if (std::find(selected_nodes.begin(), selected_nodes.end(), hap.id) != selected_nodes.end())
                    selected.emplace_back(hap);
            }
            
            std::vector<int> cooccuring_sites;
            for (auto c_mut: cooccuring)
                cooccuring_sites.emplace_back(std::stoi(c_mut.substr(0, c_mut.size()- 1)));
            
            int single_EPP_reads = 0;
            std::unordered_map<int, int> parsimony;
            for (auto rp: reads_carrying_muts) 
            {
                //Find EPPs for every read
                std::vector<haplotype> epps;
                int min_dist = std::numeric_limits<int>::max();
                std::vector<std::vector<MAT::Mutation>> epp_mut_vector; 
                for (const auto& curr_node: selected) {
                    auto curr_dist = curr_node.mutation_distance_vector(rp);
                    if ((int)curr_dist.size() <= min_dist) {
                        if ((int)curr_dist.size() < min_dist) {
                            min_dist = (int)curr_dist.size();
                            epps.clear();
                        }
                        else {
                            int curr_remaining_parsimony = 0;
                            for (auto mut: curr_dist) {
                                if (std::find(cooccuring_sites.begin(), cooccuring_sites.end(), mut.position) == cooccuring_sites.end())
                                    curr_remaining_parsimony++;
                            }

                            auto itr_epp = epps.begin();
                            while (itr_epp != epps.end()) {
                                auto prev_dist = itr_epp->mutation_distance_vector(rp);
                                int prev_remaining_parsimony = 0;
                                for (auto mut: prev_dist) {
                                    if (std::find(cooccuring_sites.begin(), cooccuring_sites.end(), mut.position) == cooccuring_sites.end())
                                        prev_remaining_parsimony++;
                                }
                                if (prev_remaining_parsimony < curr_remaining_parsimony)
                                    itr_epp = epps.erase(itr_epp);
                                else
                                    itr_epp++;
                            }   
                        }
                        epps.emplace_back(curr_node);
                        epp_mut_vector.emplace_back(curr_dist);
                    }
                }
                
                if (parsimony.find(min_dist) == parsimony.end())
                    parsimony[min_dist] = rp.degree;
                else
                    parsimony[min_dist] += rp.degree;

                if (epps.size() == 1) {
                    single_EPP_reads += rp.degree;
                    if (std::find(single_epps.begin(), single_epps.end(), epps.front().id) == single_epps.end())
                        single_epps.emplace_back(epps.front().id);
                }
                else {
                    printf("EPP:");
                    for (size_t i = 0; i < epps.size(); i++) {
                        auto node = epps[i];
                        printf(" %s(", node.id.c_str());
                        for (auto mut: epp_mut_vector[i])
                            printf("%d", mut.position);
                        printf(")"); 
                    }
                    printf("\n");
                }
            }
            printf("Single EPPs: %ld reads: %d\n", single_epps.size(), single_EPP_reads);
            for (auto epp: single_epps)
                printf("%s\n", epp.c_str());
            for (auto par: parsimony)
                printf("Parsimony: %d - %d\n", par.first, par.second);
        }
        else 
            printf("%s -> No reads cover these mutations\n", grouped_mutations.c_str());
    }
}

void arena::print_mutation_distance(const std::vector<haplotype *> &selected)
{
    std::unordered_set<int> site_read_map = this->site_read_map();
    double average_dist = 0.0;
    std::vector<std::pair<std::string, std::vector<MAT::Mutation>>>  comp_mut = this->ds.true_haplotypes();
    printf("--- ground truth %zu haps; predicted %zu haps\n", comp_mut.size(), selected.size());

    for (const auto &reference : comp_mut)
    {
        //auto sample_mutations = reference.second;
        auto sample_mutations = get_mutations(this->mat, reference.first);
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
            auto node_mutations = get_mutations(this->mat, pn->id);
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

        average_dist += (double)min_dist / comp_mut.size();
        printf("* dist: %02d true_node: %s pred_node: %s \n", min_dist, reference.first.c_str(), best_node->identifier.c_str());
    }

    printf("average mutation_distance: %0.3f\n", average_dist);
}

void arena::print_flipped_mutation_distance(const std::vector<std::pair<haplotype *, double>> &selected)
{
    std::unordered_set<int> site_read_map = this->site_read_map();
    double average_dist = 0.0;
    std::vector<std::pair<std::string, std::vector<MAT::Mutation>>>  comp_mut = this->ds.true_haplotypes();
    printf("--- ground truth %zu haps; predicted %zu haps\n", comp_mut.size(), selected.size());

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
        for (const auto &reference : comp_mut)
        {   
            //auto sample_mutations = reference.second;
            auto sample_mutations = get_mutations(this->mat, reference.first);
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
                best_node = reference.first;
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
    std::string command = "python src/WBE/sam_generation.py " + ds.directory() + " " + ds.file_prefix();
    int result = std::system(command.c_str());
    if (result)
        fprintf(stderr, "\nCannot run sam_generation.py\n");
}

void arena::resolve_unaccounted_mutations(const std::vector<std::pair<haplotype *, double>> &abundance)
{
    std::ofstream csv_haplotypes(this->ds.mutation_haplotypes_path()), csv_reads(this->ds.mutation_reads_path());
    std::string csv_print_haplotypes, csv_print_reads;

    // Read residual_mutations
    std::string file_path = this->ds.residual_mutations_path();
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
    tbb::concurrent_hash_map<std::string, std::vector<raw_read>> mutations_read_map, masked_mutations_read_map;
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
            std::vector<std::string> residual_mutations_covered, residual_mutations_masked;
            for (const auto& mut_tuple: residual_sites_covered) {
                std::string curr_residual_mut = std::to_string(std::get<0>(mut_tuple)) + std::get<1>(mut_tuple) + ":" + std::to_string(std::get<2>(mut_tuple));
                bool site_found = false;
                for (auto &mut: rp.mutations) {
                    if (mut.position == std::get<0>(mut_tuple)) {
                        // If NOT 'N' then both position and allele should match
                        if (mut.mut_nuc != 0b1111) {
                            if (MAT::get_nuc(mut.mut_nuc) == std::get<1>(mut_tuple)) {
                                // Masking the mutation on site before pushing
                                mut.mut_nuc = 0b1111;
                                residual_mutations_covered.emplace_back(curr_residual_mut);
                            }
                        }
                        // If 'N' then only position should match
                        else {
                            residual_mutations_masked.emplace_back(curr_residual_mut);
                        }
                        site_found = true;
                        break;
                    }
                }
                if ((!site_found) && (ds.reference()[std::get<0>(mut_tuple) - 1] == std::get<1>(mut_tuple))) {
                    // Add site as masked mutation
                    MAT::Mutation mut;
                    mut.ref_nuc = MAT::get_nuc_id(std::get<1>(mut_tuple));
                    mut.par_nuc = mut.ref_nuc;
                    mut.mut_nuc = 0b1111;
                    mut.position = std::get<0>(mut_tuple);
                    rp.mutations.emplace_back(mut);
                    std::sort(rp.mutations.begin(), rp.mutations.end());
                    residual_mutations_covered.emplace_back(curr_residual_mut);
                }
            }

            // Add rp to mutations_read_map if it covers residual mutations
            for (const auto& curr_residual_mut: residual_mutations_covered) {
                tbb::concurrent_hash_map<std::string, std::vector<raw_read>>::accessor ac;
                mutations_read_map.insert(ac, curr_residual_mut);
                ac->second.emplace_back(rp);
                ac.release();
            }

            // Add rp to masked_mutations_read_map if it masks residual mutations
            for (const auto& curr_residual_mut: residual_mutations_masked) {
                tbb::concurrent_hash_map<std::string, std::vector<raw_read>>::accessor ac;
                masked_mutations_read_map.insert(ac, curr_residual_mut);
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

    // Add masked_mutations_read_map to mutations_read_map for residual mutation to haplotype mapping
    for (const auto &m_r : masked_mutations_read_map)
    {
        tbb::concurrent_hash_map<std::string, std::vector<raw_read>>::accessor ac;
        mutations_read_map.insert(ac, m_r.first);
        ac->second.insert(ac->second.end(), m_r.second.begin(), m_r.second.end());
        ac.release();
    }
    masked_mutations_read_map.clear();

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
                    int min_dist = std::numeric_limits<int>::max();
                    for (const auto& curr_node_abun: abundance) {
                        int curr_dist = curr_node_abun.first->mutation_distance(rp);
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

void arena::dump_haplotypes(const std::vector<std::pair<haplotype *, double>> &abundance)
{
    std::ofstream tsv(this->ds.haplotype_tsv_path());
    std::string tsv_print;
    
    for (size_t i = 0; i < abundance.size(); i++)
    {
        auto hap = abundance[i].first;
        auto hap_mutations = ::get_mutations(mat, hap->id);
        std::sort(hap_mutations.begin(), hap_mutations.end());
        int start_idx = 1;
        std::string md = "MD:Z:";
        auto hap_sequence = this->reference();
        for (const auto& mut: hap_mutations)
        {
            hap_sequence[mut.position - 1] = MAT::get_nuc(mut.mut_nuc);
            md += std::to_string(mut.position - start_idx) + MAT::get_nuc(mut.ref_nuc);
            start_idx = mut.position + 1;
        }
        if (hap_sequence.size() - start_idx + 1)
            md += std::to_string(hap_sequence.size() - start_idx + 1);
        tsv_print = hap->id + "\t0\tNC_045512.2\t1\t60\t" + std::to_string(this->reference().size()) + "M\t*\t0\t0\t" + hap_sequence + "\t" + std::string(hap_sequence.length(), '?') + "\t" + md + "\tRG:Z:group\n";
        tsv << tsv_print;
    }
}