#ifndef REDSET_INTERNAL_H
#define REDSET_INTERNAL_H

#include "kvtree.h"

#include "redset.h"
#include "redset_io.h"
#include "redset_util.h"
#include "redset_lofi.h"

#define REDSET_VERSION "1.0"

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
  int       replicas;       /* number of partner replicas */
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
  kvtree*   group_map;       /* kvtree that maps group rank to world rank */
  int       lhs_rank;        /* rank which is one less (with wrap to highest) within set */
  int       lhs_rank_world;  /* rank of lhs process in comm world */
  char*     lhs_hostname;    /* hostname of lhs process */
  int       rhs_rank;        /* rank which is one more (with wrap to lowest) within set */
  int       rhs_rank_world;  /* rank of rhs process in comm world */
  char*     rhs_hostname;    /* hostname of rhs process */
  int       encoding;        /* number of encoding blocks */
  int       gf_bits;         /* number of bits in Galois Field */
  int       gf_size;         /* number of elements in Galois Field */
  unsigned int*  gf_log;     /* log_2 table */
  unsigned int*  gf_exp;     /* exp_2 table (inverse log) */
  unsigned int*  gf_imult;   /* computes multiplicative inverse for each value */
  unsigned char* gf_premult; /* pre-computed products for all elements against a constant */
  unsigned int*  mat;        /* encoding matrix (ranks + encoding) x ranks */
} redset_reedsolomon;

typedef struct {
  int count;
  const char** files;
} redset_list;

int redset_set_partners(
  MPI_Comm parent_comm, MPI_Comm comm, int dist,
  int* lhs_rank, int* lhs_rank_world, char** lhs_hostname,
  int* rhs_rank, int* rhs_rank_world, char** rhs_hostname
);

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

int redset_construct_partner(
  MPI_Comm parent_comm,
  redset_base* d,
  int replicas
);

int redset_construct_xor(
  MPI_Comm parent_comm,
  redset_base* d
);

int redset_construct_rs(
  MPI_Comm parent_comm,
  redset_base* d,
  int encoding
);

int redset_delete_partner(
  redset_base* d
);

int redset_delete_xor(
  redset_base* d
);

int redset_delete_rs(
  redset_base* d
);

int redset_store_to_kvtree_partner(
  const redset_base* d,
  kvtree* hash
);

int redset_read_from_kvtree_partner(
  const kvtree* hash,
  int* outreplicas
);

int redset_store_to_kvtree_rs(
  const redset_base* d,
  kvtree* hash
);

int redset_read_from_kvtree_rs(
  const kvtree* hash,
  int* outencoding
);

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

int redset_encode_reddesc_rs(
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

int redset_apply_rs(
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

int redset_recover_rs(
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

int redset_unapply_rs(
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

redset_list* redset_filelist_get_rs(
  const char* name,
  redset_base* d
);

#endif /* REDSET_INTERNAL_H */
