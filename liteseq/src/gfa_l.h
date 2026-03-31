#ifndef LQ_GFA_L_H
#define LQ_GFA_L_H

#include "../include/liteseq/gfa.h"

struct l_thread_meta {
	edge *edges;
	line *l_lines;
	idx_t l_line_count;
};

void *t_handle_l(void *l_meta);
#endif // LQ_GFA_L_H
