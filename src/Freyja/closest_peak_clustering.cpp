#include "src/usher_graph.hpp"

// Define a structure to hold abs_value (float) and a boolean vector
struct costMutations {
    double cost;                
    std::vector<int> mutations; 
};

// Function to read CSV file and return a vector of costMutations structs
std::vector<struct costMutations> readCSV(const std::string& filename) {
    std::vector<struct costMutations> cost_mutations;
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Could not open the file " << filename << std::endl;
        return cost_mutations;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string token;
        struct costMutations entry;

        // Read the first value (abs_value) as a float
        std::getline(ss, token, ',');
        entry.cost = std::stod(token);

        // Read the rest of the tokens as boolean values
        while (std::getline(ss, token, ',')) {
            entry.mutations.emplace_back(std::stoi(token));
        }

        // Add the entry to the vector of costMutations structs
        cost_mutations.emplace_back(entry);
    }

    file.close();
    return cost_mutations;
}

int main(int argc, char* argv[]) {
    // Check if the command line argument for epsilon is provided
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <epsilon>" << std::endl;
        return 1;
    }
    // Get epsilon from the command line and convert it to float
    float eps = std::stof(argv[1]);

    // Read the CSV file and get the vector of costMutations structs
    std::string input_filename = "./closest_peak_search.csv";
    std::string output_filename = "./peaks_clustered.csv";
    std::vector<struct costMutations> data = readCSV(input_filename);
    if (data.empty())
        return 1;

    // Create indices vector and sort it
    std::vector<int> indices(data.size()), closest_indices(data.size());
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
            auto curr_mutations = curr_data.mutations;
            int min_dist = INT_MAX;
            int closest_j = INT_MAX;
            for (size_t j = i+1; j < indices.size(); j++) 
            {
                int jdx = indices[j];
                auto cmp_data = data[jdx];
                auto cmp_mutations = cmp_data.mutations;
                int curr_dist = 0;
                for (size_t k = 0; k < curr_mutations.size(); k++) {
                    if (curr_mutations[k] != cmp_mutations[k])
                        curr_dist++;
                }
                if (curr_dist < min_dist) {
                    min_dist = curr_dist;
                    closest_j = jdx;
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
            int jdx = closest_indices[i];
            data[jdx].cost += data[idx].cost;
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
