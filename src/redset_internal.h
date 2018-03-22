#ifndef REDSET_INTERNAL_H
#define REDSET_INTERNAL_H

#include "kvtree.h"
#include "redset.h"

int redset_encode_reddesc_single(
  kvtree* hash,
  const char* name,
  const redset* d
);

int redset_encode_reddesc_partner(
  kvtree* hash,
  const char* name,
  const redset* d
);

int redset_encode_reddesc_xor(
  kvtree* hash,
  const char* name,
  const redset* d
);


int redset_apply_single(
  int numfiles,
  const char** files,
  const char* name,
  const redset* d
);

int redset_apply_partner(
  int numfiles,
  const char** files,
  const char* name,
  const redset* d
);

int redset_apply_xor(
  int numfiles,
  const char** files,
  const char* name,
  const redset* d
);


int redset_recover_single(
  const char* name,
  const redset* d
);

int redset_recover_partner(
  const char* name,
  const redset* d
);

int redset_recover_xor(
  const char* name,
  const redset* d
);


int redset_unapply_single(
  const char* name,
  const redset* d
);

int redset_unapply_partner(
  const char* name,
  const redset* d
);

int redset_unapply_xor(
  const char* name,
  const redset* d
);


redset_filelist* redset_filelist_get_single(
  const char* name,
  redset* d
);

redset_filelist* redset_filelist_get_partner(
  const char* name,
  redset* d
);

redset_filelist* redset_filelist_get_xor(
  const char* name,
  redset* d
);

#endif /* REDSET_INTERNAL_H */
