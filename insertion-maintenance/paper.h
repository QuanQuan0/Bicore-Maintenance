#pragma once
#ifndef PAPER_H
#define PAPER_H
#include "bigraph.h"
#include "dyn_rebuild.h"
#include <fstream>

extern bool DELETEC;
extern bool APC;
extern long long DELETECOUNT;
extern long long APCOUNT;

extern void alphaCopyPeel_for_kcore(int left_k, BiGraph& g);

extern void __coreIndexBasic(BiGraph& g);

extern void coreIndexBasic(BiGraph& g);

extern int coreIndexKCore(BiGraph& g);

extern int coreIndexKCore_with_pruning_rules(BiGraph& g);

extern void findabcore_core_index(int alpha, int beta, BiGraph& g, std::vector<bool>& left_sub, std::vector<bool>& right_sub);

extern bool test_abcore_equal(std::vector<bool>& left1, std::vector<bool>& right1, std::vector<bool>& left2, std::vector<bool>& right2);

extern int maxkcorenumber(BiGraph& g);

#endif // !GEPHI_H