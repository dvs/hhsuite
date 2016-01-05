#ifndef _DOTHELIX_H_

#include <vector>

typedef struct _diagonal_segment {
    int start;
    int end;
    double value;
} diagonal_segment;

extern int check_by_full_search;

void find_all_segments(std::vector<diagonal_segment> * r, std::vector<double> v, int start, int end, double M, double D, double threshold);
void sort_segments(std::vector<diagonal_segment> * r);

#endif
