#ifndef LQ_GFA_S_H
#define LQ_GFA_S_H

#include "../include/liteseq/gfa.h"

struct s_thread_meta {
	vtx **vertices;
	line *s_lines;
	idx_t s_line_count;
	bool inc_vtx_labels;
};

void *t_handle_s(void *s_meta);
#endif // LQ_GFA_S_H
