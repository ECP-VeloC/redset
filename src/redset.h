/*
 * Copyright (c) 2009, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * Written by Adam Moody <moody20@llnl.gov>.
 * LLNL-CODE-411039.
 * All rights reserved.
 * This file is part of The Scalable Checkpoint / Restart (SCR) library.
 * For details, see https://sourceforge.net/projects/scalablecr/
 * Please also read this file: LICENSE.TXT.
*/

#ifndef REDSET_H
#define REDSET_H

#include "mpi.h"
#include "kvtree.h"

#define REDSET_SUCCESS (0)

#define REDSET_COPY_NULL    (0)
#define REDSET_COPY_SINGLE  (1)
#define REDSET_COPY_PARTNER (2)
#define REDSET_COPY_XOR     (3)

/*
=========================================
Define redundancy descriptor structure
=========================================
*/

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
} redset;

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
} redset_filelist;

/*
=========================================
Redundancy descriptor functions
=========================================
*/

/* initialize library */
int redset_init();

/* shutdown library */
int redset_finalize();

/* create a new redundancy set descriptor */
int redset_create(
  int type,          /* redundancy encoding type */
  MPI_Comm comm,     /* process group participating in set */
  const char* group, /* string specifying procs in the same failure group */
  redset* d          /* output redundancy descriptor */
);

/* free any memory associated with the specified redundancy descriptor */
int redset_delete(
  redset* d
);

/* apply redundancy scheme to file and return number of bytes copied
 * in bytes parameter */
int redset_apply(
  int numfiles,
  const char** files,
  const char* name,
  const redset* d
);

/* rebuilds files for specified dataset id using specified redundancy descriptor,
 * adds them to filemap, and returns REDSET_SUCCESS if all processes succeeded */
int redset_recover(
  MPI_Comm comm,
  const char* name,
  redset* d
);

/* deletes redundancy data that was added in redset_apply,
 * which is useful when cleaning up */
int redset_unapply(
  const char* name,
  redset* d
);

/* return list of files added by redundancy scheme */
redset_filelist* redset_filelist_get(
  const char* name,
  redset* d
);

/* free file list allocated by call to redset_filelist */
int redset_filelist_release(redset_filelist** plist);

/* return number of files in file list */
int redset_filelist_count(redset_filelist* list);

/* return name of file at specified index, where
 * 0 <= index < count and count is value returned by redset_filelist_count */
const char* redset_filelist_file(redset_filelist* list, int index);

#define ER_SUCCESS (0)
#define ER_FAILURE (1)

#define ER_DIRECTION_ENCODE  (1)
#define ER_DIRECTION_REBUILD (2)

int ER_Init(const char* conf_file);

int ER_Finalize();

int ER_Create_Scheme(
  MPI_Comm comm,
  const char* failure_domain,
  int encoding_blocks,
  int erasure_blocks
);

int ER_Free_Scheme(int scheme_id);

/* create a named set, and specify whether it should be encoded or recovered */
int ER_Create(
  const char* name,
  int direction,
  int scheme_id /* encoding scheme to be applied */
);

/* adds file to specified set id */
int ER_Add(
  int set_id,
  const char* file
);

/* initiate encode/rebuild operation on specified set id */
int ER_Dispatch(
  int set_id
);

/* tests whether ongoing dispatch operation to finish,
 * returns 1 if done, 0 otherwise */
int ER_Test(
  int set_id
);

/* wait for ongoing dispatch operation to finish */
int ER_Wait(
  int set_id
);

/* free internal resources associated with set id */
int ER_Free(
  int set_id
);

#endif /* REDSET_H */
