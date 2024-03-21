#include "pwalign.h"

#define CALLMETHOD_DEF(fun, numArgs) {#fun, (DL_FUNC) &fun, numArgs}

static const R_CallMethodDef callMethods[] = {

/* align_utils.c */
	CALLMETHOD_DEF(PairwiseAlignments_nmatch, 4),
	CALLMETHOD_DEF(AlignedXStringSet_nchar, 1),
	CALLMETHOD_DEF(AlignedXStringSet_align_aligned, 2),
	CALLMETHOD_DEF(PairwiseAlignmentsSingleSubject_align_aligned, 3),
	CALLMETHOD_DEF(align_compareStrings, 6),

/* align_pairwiseAlignment.c */
	CALLMETHOD_DEF(XStringSet_align_pairwiseAlignment, 14),
	CALLMETHOD_DEF(XStringSet_align_distance, 12),

	{NULL, NULL, 0}
};

void R_init_pwalign(DllInfo *info)
{
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
	R_useDynamicSymbols(info, 0);
	return;
}

