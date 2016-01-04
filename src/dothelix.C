#include <cmath>
#include <algorithm>
#include <iostream>

#include "dothelix.h"

int check_by_full_search = 0;

static inline double measure(double sum, int l, double M, double D)
{
    return (sum - (M + 0.0) * l) / (D * sqrt(l));
}

static diagonal_segment find_max_segment_reverse(std::vector<double> v, int start0, int end0, double M, double D);

static diagonal_segment find_max_segment_full_search(std::vector<double> v, int start, int end, double M, double D)
{
    int L = end - start + 1;
    double sums[L];
    double sum = 0.0;
    sums[0] = 0.0;
    for (int i = start; i <= end; i++) {
        sum += v[i];
        sums[i - start + 1] = sum;
    }
    double max_value = 0.0;
    int max_start = -1;
    int max_end = -1;
    for(int s = start; s <= end; s++) {
        for(int e = s; e <= end; e++) {
            int l = e - s + 1;
            double sum = sums[e - start + 1] - sums[s - start];
            double value = measure(sum, l, M, D);
            if (max_value < value) {
                max_value = value;
                max_start = s;
                max_end = e;
            }
        }
    }
    diagonal_segment result;
    result.start = max_start;
    result.end = max_end;
    result.value = max_value;
    return result;
}

static diagonal_segment find_max_segment(std::vector<double> v, int start0, int end0, double M, double D)
{
    int start = start0;
    int end = end0;
    diagonal_segment result;
    result.start = result.end = -1;
    result.value = 0.0;

    //! dbg
    // printf("start = %d, end = %d\n", start, end);
    if (start > end) {
        return result;
    }
    while (v[start] < M && start <= end) {
        start++;
    }
    if (start == end) {
        double value = (v[start] - M) / D;
        if (value > 0) {
            result.start = result.end = start;
            result.value = value;
        }
        return result;
    }

    int first_negative_index = -1;
    double sum = 0.0;
    int max_positive_index = -1;
    double max_positive_value = 0.0;

    for (int i = start; i <= end; i++) {
        sum += v[i];
        int l = i - start + 1;
        double value = measure(sum, l, M, D);
        if (value < 0 && i > start && first_negative_index == -1) {
            first_negative_index = i;
        }
        if (max_positive_value < value) {
            max_positive_value = value;
            max_positive_index = i;
        }
    }
    if (max_positive_index == -1) {
        return result;
    }

    int split = 1;
    int start1, end1, start2, end2;
    if (first_negative_index >= 0) { // segment is split by first_negative_index
        start1 = start;
        end1 = first_negative_index - 1;
        start2 = first_negative_index;
        end2 = end;
    } else {
        if (max_positive_index == end) { // we cannot split
            split = 0;
        } else { // segment is split by max_positive_index
            start1 = start;
            end1 = max_positive_index;
            start2 = max_positive_index + 1;
            end2 = end;
        }
    }

    if (split) {
        //! dbg
        // printf("first_negative_index = %d, max_positive_index = %d\n", first_negative_index, max_positive_index);
        diagonal_segment result1 = find_max_segment(v, start1, end1, M, D);
        diagonal_segment result2 = find_max_segment(v, start2, end2, M, D);
        if (result1.value > result2.value) {
            return result1;
        } else {
            return result2;
        }
    } else {
        result = find_max_segment_reverse(v, start, end, M, D);
        return result;
    }
}

static diagonal_segment find_max_segment_reverse(std::vector<double> v, int start0, int end0, double M, double D)
{
    int start = start0;
    int end = end0;
    diagonal_segment result;
    result.start = result.end = -1;
    result.value = 0.0;

    //! dbg
    // printf("REVERSE: start = %d, end = %d\n", start, end);
    if (start > end) {
        return result;
    }
    while (v[end] < M && start <= end) {
        end--;
    }
    if (start == end) {
        double value = (v[start] - M) / D;
        if (value > 0) {
            result.start = result.end = start;
            result.value = value;
        }
        return result;
    }

    int first_negative_index = -1;
    double sum = 0.0;
    int max_positive_index = -1;
    double max_positive_value = 0.0;

    for (int i = end; start <= i; i--) {
        sum += v[i];
        int l = end - i + 1;
        double value = measure(sum, l, M, D);
        if (value < 0 && i < end && first_negative_index == -1) {
            first_negative_index = i;
        }
        if (max_positive_value < value) {
            max_positive_value = value;
            max_positive_index = i;
        }
    }
    if (max_positive_index == -1) {
        return result;
    }

    int split = 1;
    int start1, end1, start2, end2;
    if (first_negative_index >= 0) { // segment is split by first_negative_index
        start1 = start;
        end1 = first_negative_index;
        start2 = first_negative_index + 1;
        end2 = end;
    } else {
        if (max_positive_index == start) { // we cannot split
            split = 0;
        } else { // segment is split by max_positive_index
            start1 = start;
            end1 = max_positive_index - 1;
            start2 = max_positive_index;
            end2 = end;
        }
    }

    if (split) {
        //! dbg
        // printf("first_negative_index = %d, max_positive_index = %d\n", first_negative_index, max_positive_index);
        diagonal_segment result1 = find_max_segment_reverse(v, start1, end1, M, D);
        diagonal_segment result2 = find_max_segment_reverse(v, start2, end2, M, D);
        if (result1.value > result2.value) {
            return result1;
        } else {
            return result2;
        }
    } else {
        if (start != start0 || end != end0) {
            result = find_max_segment(v, start, end, M, D);
        } else {
            result = find_max_segment_full_search(v, start, end, M, D);
        }
        return result;
    }
}


bool _compare_segments(diagonal_segment a, diagonal_segment b) {return a.value > b.value;}

void sort_segments(std::vector<diagonal_segment> * a)
{
    std::sort(a->begin(), a->end(), _compare_segments);
}

void find_all_segments(std::vector<diagonal_segment> * r, std::vector<double> v, int start, int end, double M, double D, double threshold)
{
    diagonal_segment s = find_max_segment(v, start, end, M, D);
    if (check_by_full_search) {
        diagonal_segment s0 = find_max_segment_full_search(v, start, end, M, D);
        if (s.start != s0.start || s.end != s0.end) {
            std::cerr << "Wrong segment from find_max_segment(): [" << s.start << "; " << s.end << "]" << std::endl;
            s = s0;
        }
    }
    if (s.start == -1) {
        return;
    }
    if (s.value >= threshold) {
        r->push_back(s);
    }
    if (s.start > start) {
        find_all_segments(r, v, start, s.start - 1, M, D, threshold);
    }
    if (s.end < end) {
        find_all_segments(r, v, s.end + 1, end, M, D, threshold);
    }
}


