#include <vector>
#include <tbb/blocked_range.h>
#include <tbb/queuing_mutex.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_sort.h>
#include <iostream>
#include <fstream>
#include <sstream>

// Define a structure to hold abs_value (float) and a boolean vector
struct costMutations {
    double cost;                
    std::vector<int> mutation_idx; 
};

struct Mutation {
    int pos;
    char mut;
};

// Function to read CSV file and return a vector of costMutations structs
std::vector<costMutations> readCSV(const std::string& filename, std::vector<Mutation>& mutation_list) {
    std::vector<costMutations> cost_mutations;
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Could not open the file " << filename << std::endl;
        return cost_mutations;
    }

    std::string line;
    // Read the header line
    if (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string token;

        // Skip the first empty column
        std::getline(ss, token, ',');

        // Read the mutation headers
        while (std::getline(ss, token, ',')) {
            Mutation mut;
            mut.pos = std::stoi(token.substr(1, token.size() - 2));
            mut.mut = token.back();
            mutation_list.emplace_back(mut);
        }
    }

    // Process the remaining lines
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string token;
        costMutations entry;

        // Read the first value (abs_value) as a float
        std::getline(ss, token, ',');
        entry.cost = std::stod(token);

        // Read the rest of the tokens as boolean values
        int idx = 0; 
        while (std::getline(ss, token, ',')) {
            if (std::stoi(token) != 0) { 
                entry.mutation_idx.emplace_back(idx);
            }
            idx++;
        }

        // Add the entry to the vector of costMutations structs
        cost_mutations.emplace_back(entry);
    }

    file.close();
    return cost_mutations;
}

int mutation_distance(const std::vector<Mutation>& node1_mutations, const std::vector<Mutation>& node2_mutations) {
    int muts = 0, i = 0, j = 0;

    while (i < (int) node1_mutations.size() || j < (int) node2_mutations.size()) {
        if (i == (int) node1_mutations.size()) {                
            ++muts;
            ++j;
        }
        else if (j == (int) node2_mutations.size()) {
            ++muts;
            ++i;
        }
        else if (node1_mutations[i].pos < node2_mutations[j].pos) {
            ++muts;
            ++i;
        }
        else if (node1_mutations[i].pos > node2_mutations[j].pos) {
            ++muts;
            ++j;
        }
        else if (node1_mutations[i].pos == node2_mutations[j].pos && node1_mutations[i].mut != node2_mutations[j].mut) {
            ++muts;
            ++i; ++j;
        }
        else {
            ++i; ++j;
        }
    }
    
    return muts;
}

int main(int argc, char* argv[]) {
    // Check if the command line argument for epsilon is provided
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <epsilon> <file_path>" << std::endl;
        return 1;
    }
    // Get epsilon from the command line and convert it to float
    float eps = std::stof(argv[1]);
    std::string file_path = argv[2];

    // Read the CSV file and get the vector of costMutations structs
    std::string input_filename = file_path + "/closest_peak_search.csv";
    std::string output_filename = file_path + "/peaks_clustered.csv";
    std::vector<Mutation> mutation_list;
    std::vector<costMutations> data = readCSV(input_filename, mutation_list);
    if (data.empty())
        return 1;

    // Create indices vector and sort it
    std::vector<int> indices(data.size());
    std::vector<std::vector<int>> closest_indices(data.size());
    for (size_t i = 0; i < data.size(); ++i) {
        indices[i] = i;
    }
    tbb::parallel_sort(indices.begin(), indices.end(), [&data](const int& i1, const int& i2) {
        return data[i1].cost < data[i2].cost;  
    });

    // Modify the cost values
    tbb::queuing_mutex my_mutex;
    static tbb::affinity_partitioner ap;
    tbb::parallel_for(tbb::blocked_range<size_t>(0, indices.size()), [&](const tbb::blocked_range<size_t>& r) 
    {
        for (size_t i = r.begin(); i < r.end(); ++i) 
        { 
            int idx = indices[i];
            auto curr_data = data[idx];
            std::vector<Mutation> curr_mutations;
            for (int idx: curr_data.mutation_idx)
                curr_mutations.emplace_back(mutation_list[idx]);
            int min_dist = INT_MAX;
            std::vector<int> closest_j;
            for (size_t j = i+1; j < indices.size(); j++) 
            {
                int jdx = indices[j];
                auto cmp_data = data[jdx];
                std::vector<Mutation> cmp_mutations;
                for (int idx: cmp_data.mutation_idx)
                    cmp_mutations.emplace_back(mutation_list[idx]);
                int curr_dist = mutation_distance(curr_mutations, cmp_mutations);
                if (curr_dist < min_dist) {
                    min_dist = curr_dist;
                    closest_j.clear();
                    closest_j.emplace_back(jdx);
                }
                else if (curr_dist == min_dist) {
                    closest_j.emplace_back(jdx);
                }
            }
            {
                tbb::queuing_mutex::scoped_lock my_lock{my_mutex};
                closest_indices[i] = closest_j;
            }
        }
    });

    // Updating cost values
    for (size_t i = 0; i < indices.size(); ++i) 
    {
        int idx = indices[i];
        if ((data[idx].cost < eps) && (i+1 < indices.size())) {
            double inc_cost = data[idx].cost / (double) closest_indices[i].size();
            for (int jdx: closest_indices[i]) {
                data[jdx].cost += inc_cost;
            }
            data[idx].cost = 0;
        }
    }
                
    // Open File for writing
    std::ofstream outputFile(output_filename);
    if (!outputFile.is_open()) {
        std::cerr << "Error: Could not open the file for writing!" << std::endl;
        return 1;
    }

    // Write the cost in the file based on original order of data
    for (size_t i = 0; i < data.size(); ++i)
        outputFile << data[i].cost << std::endl;
    outputFile.close();
    
    return 0;
}