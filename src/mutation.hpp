#pragma once

#include <cstring>
#include <cctype>

enum nuc {
    NUC_A   = 0b00001,
    NUC_C   = 0b00010,
    NUC_G   = 0b00100,
    NUC_T   = 0b01000,
    NUC_GAP = 0b10000,
    NUC_N   = 0b11111,
};

// does not handle gap case
static nuc nuc_from_pannuc(int pan_nuc) {
    switch (pan_nuc) {
    case 1:
        return NUC_A;
    case 2:
        return NUC_C;
    case 4:
        return NUC_G;
    case 8:
        return NUC_T;
    default:
        printf("Unknown nucleotide conversion pan_nuc %d\n", pan_nuc);
        return NUC_N;
    }
}

static nuc nuc_from_char(char c) {
    switch (tolower(c)) {
    case 'a':
        return NUC_A;
    case 'c':
        return NUC_C;
    case 'g':
        return NUC_G;
    case 't':
        return NUC_T;
    case '_':
        return NUC_GAP;
    case 'n':
    default:
        return NUC_N;
    }
}

static char char_from_nuc(nuc n) {
    switch (n) {
    case NUC_A:
        return 'A';
    case NUC_C:
        return 'C';
    case NUC_G:
        return 'G';
    case NUC_T:
        return 'T';
    case NUC_GAP:
        return '_';
    case NUC_N:
    default:
        return 'N';
    }
}

struct mutation {
    int pos;
    nuc ref, mut;

    // note that it says if it's an indel relative to reference
    // not relative to parent
    bool is_indel() const {
        return this->ref == NUC_GAP || this->mut == NUC_GAP;
    }

    std::string get_string() const {
        return char_from_nuc(ref) + std::to_string(pos) + char_from_nuc(mut);
    }

    bool operator< (const mutation& other) const {
        return this->pos < other.pos;
    }
};