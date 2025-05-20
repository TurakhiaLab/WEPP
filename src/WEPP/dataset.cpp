#include "dataset.hpp"
#include "sam2pb.hpp"

std::vector<raw_read> dataset::reads() const {
    std::string protoname = this->pb_path();
    std::unordered_map<std::string, std::vector<std::string>> reverse;
    return load_reads_from_proto(this->reference(), protoname, reverse);
}

std::unordered_map<std::string, std::vector<std::string>> dataset::read_reverse_merge() const {
    std::string protoname = this->pb_path();
    std::unordered_map<std::string, std::vector<std::string>> reverse;
    load_reads_from_proto(this->reference(), protoname, reverse);
    return reverse;
}