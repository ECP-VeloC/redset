#ifndef REDSET_INTERNAL_H
#define REDSET_INTERNAL_H

#include "kvtree.h"
#include "redset.h"

typedef struct {
  int      enabled;        /* flag indicating whether this descriptor is active */
  int      type;           /* redundancy scheme to apply */
  void*    state;          /* pointer to extra state depending on copy type */
  MPI_Comm parent_comm;    /* parent communicator */
  MPI_Comm comm;           /* communicator holding procs for this scheme */
  int      groups;         /* number of redundancy sets */
  int      group_id;       /* unique id assigned to this redundancy set */
  int      ranks;          /* number of ranks in this set */
  int      rank;           /* caller's rank within its set */
} redset_base;

typedef struct {
  int       lhs_rank;       /* rank which is one less (with wrap to highest) within set */
  int       lhs_rank_world; /* rank of lhs process in comm world */
  char*     lhs_hostname;   /* hostname of lhs process */
  int       rhs_rank;       /* rank which is one more (with wrap to lowest) within set */
  int       rhs_rank_world; /* rank of rhs process in comm world */
  char*     rhs_hostname;   /* hostname of rhs process */
} redset_partner;

typedef struct {
  kvtree*   group_map;      /* kvtree that maps group rank to world rank */
  int       lhs_rank;       /* rank which is one less (with wrap to highest) within set */
  int       lhs_rank_world; /* rank of lhs process in comm world */
  char*     lhs_hostname;   /* hostname of lhs process */
  int       rhs_rank;       /* rank which is one more (with wrap to lowest) within set */
  int       rhs_rank_world; /* rank of rhs process in comm world */
  char*     rhs_hostname;   /* hostname of rhs process */
} redset_xor;

typedef struct {
  int count;
  const char** files;
} redset_list;

/* convert the specified redundancy descritpor into a corresponding
 * kvtree */
int redset_store_to_kvtree(
  const redset_base* d,
  kvtree* kv
);

/* build a redundancy descriptor corresponding to the specified kvtree,
 * this function is collective, it differs from create_from_kvtree in
 * that it uses group id and group rank values to restore a descriptor
 * that was previously created */
int redset_restore_from_kvtree(
  const kvtree* kv,
  redset_base* d
);

/* capture file metadata for file into meta */
int redset_meta_encode(const char* file, kvtree* meta);

/* apply file metadata in meta to file */
int redset_meta_apply(const char* file, const kvtree* meta);

int redset_encode_reddesc_single(
  kvtree* hash,
  const char* name,
  const redset_base* d
);

int redset_encode_reddesc_partner(
  kvtree* hash,
  const char* name,
  const redset_base* d
);

int redset_encode_reddesc_xor(
  kvtree* hash,
  const char* name,
  const redset_base* d
);


int redset_apply_single(
  int numfiles,
  const char** files,
  const char* name,
  const redset_base* d
);

int redset_apply_partner(
  int numfiles,
  const char** files,
  const char* name,
  const redset_base* d
);

int redset_apply_xor(
  int numfiles,
  const char** files,
  const char* name,
  const redset_base* d
);


int redset_recover_single(
  const char* name,
  const redset_base* d
);

int redset_recover_partner(
  const char* name,
  const redset_base* d
);

int redset_recover_xor(
  const char* name,
  const redset_base* d
);


int redset_unapply_single(
  const char* name,
  const redset_base* d
);

int redset_unapply_partner(
  const char* name,
  const redset_base* d
);

int redset_unapply_xor(
  const char* name,
  const redset_base* d
);


redset_list* redset_filelist_get_single(
  const char* name,
  redset_base* d
);

redset_list* redset_filelist_get_partner(
  const char* name,
  redset_base* d
);

redset_list* redset_filelist_get_xor(
  const char* name,
  redset_base* d
);

#endif /* REDSET_INTERNAL_H */
