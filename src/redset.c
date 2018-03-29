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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <limits.h>
#include <unistd.h>

#include "mpi.h"
#include "mpi_rank_str.h"

#include "kvtree.h"
#include "kvtree_util.h"

#include "redset_util.h"
#include "redset.h"
#include "redset_internal.h"

int redset_init()
{
  /* read our hostname */
  char hostname[HOST_NAME_MAX + 1];
  gethostname(hostname, sizeof(hostname));
  redset_hostname = strdup(hostname);

  /* get page size */
  redset_page_size = sysconf(_SC_PAGESIZE);

  /* set MPI buffer size */
  redset_mpi_buf_size = 1024 * 1024;

  /* TODO: allow caller to specify process group some other way */
  /* duplicate communicator */
  MPI_Comm_dup(MPI_COMM_WORLD, &redset_comm);

  /* set our global rank */
  MPI_Comm_rank(redset_comm, &redset_rank);

  return REDSET_SUCCESS;
}

int redset_finalize()
{
  redset_free(&redset_hostname);
  MPI_Comm_free(&redset_comm);
  return REDSET_SUCCESS;
}

/* given a comm as input, find the left and right partner
 * ranks and hostnames */
int redset_set_partners(
  MPI_Comm parent_comm, MPI_Comm comm, int dist,
  int* lhs_rank, int* lhs_rank_world, char** lhs_hostname,
  int* rhs_rank, int* rhs_rank_world, char** rhs_hostname)
{
  /* find our position in the communicator */
  int my_rank, ranks;
  MPI_Comm_rank(comm, &my_rank);
  MPI_Comm_size(comm, &ranks);

  /* shift parter distance to a valid range */
  while (dist > ranks) {
    dist -= ranks;
  }
  while (dist < 0) {
    dist += ranks;
  }

  /* compute ranks to our left and right partners */
  int lhs = (my_rank + ranks - dist) % ranks;
  int rhs = (my_rank + ranks + dist) % ranks;
  (*lhs_rank) = lhs;
  (*rhs_rank) = rhs;

  /* shift hostnames to the right */
  redset_str_sendrecv(redset_hostname, rhs, lhs_hostname, lhs, comm);

  /* shift hostnames to the left */
  redset_str_sendrecv(redset_hostname, lhs, rhs_hostname, rhs, comm);

  MPI_Request request[2];
  MPI_Status  status[2];

  /* get our rank in parent communicator */
  int parent_rank;
  MPI_Comm_rank(parent_comm, &parent_rank);

  /* shift rank in comm_world to the right */
  MPI_Irecv(lhs_rank_world, 1, MPI_INT, lhs, 0, comm, &request[0]);
  MPI_Isend(&parent_rank,    1, MPI_INT, rhs, 0, comm, &request[1]);
  MPI_Waitall(2, request, status);

  /* shift rank in comm_world to the left */
  MPI_Irecv(rhs_rank_world, 1, MPI_INT, rhs, 0, comm, &request[0]);
  MPI_Isend(&parent_rank,    1, MPI_INT, lhs, 0, comm, &request[1]);
  MPI_Waitall(2, request, status);

  return REDSET_SUCCESS;
}

/*
=========================================
Redundancy descriptor functions
=========================================
*/

/* initialize the specified redundancy descriptor */
static int redset_initialize(redset* d)
{
  /* check that we got a valid redundancy descriptor */
  if (d == NULL) {
    printf("No redundancy descriptor to fill from hash @ %s:%d",
      __FILE__, __LINE__
    );
    return REDSET_FAILURE;
  }

  /* initialize the descriptor */
  d->enabled        =  0;
  d->type           = REDSET_COPY_NULL;
  d->state          = NULL;
  d->comm           = MPI_COMM_NULL;
  d->groups         =  0;
  d->group_id       = -1;
  d->ranks          =  0;
  d->rank           = MPI_PROC_NULL;

  return REDSET_SUCCESS;
}

/* given a redundancy descriptor with all top level fields filled in
 * allocate and fill in structure for partner specific fields in state */
static int redset_create_partner(MPI_Comm parent_comm, redset* d)
{
  int rc = REDSET_SUCCESS;

  /* allocate a new structure to hold partner state */
  redset_partner* state = (redset_partner*) REDSET_MALLOC(sizeof(redset_partner));

  /* attach structure to reddesc */
  d->state = (void*) state;

  /* record group rank, world rank, and hostname of left and right partners */
  redset_set_partners(
    parent_comm, d->comm, 1,
    &state->lhs_rank, &state->lhs_rank_world, &state->lhs_hostname,
    &state->rhs_rank, &state->rhs_rank_world, &state->rhs_hostname
  );

  /* check that we got valid partners */
  if (state->lhs_hostname == NULL ||
      state->rhs_hostname == NULL ||
      strcmp(state->lhs_hostname, "") == 0 ||
      strcmp(state->rhs_hostname, "") == 0 ||
      strcmp(state->lhs_hostname, redset_hostname) == 0 ||
      strcmp(state->rhs_hostname, redset_hostname) == 0)
  {
    /* disable this descriptor */
    d->enabled = 0;
    redset_warn("Failed to find partner processes for redundancy descriptor, disabling, too few nodes? @ %s:%d",
      __FILE__, __LINE__
    );
    rc = REDSET_FAILURE;
  } else {
    redset_dbg(2, "LHS partner: %s (%d)  -->  My name: %s (%d)  -->  RHS partner: %s (%d)",
      state->lhs_hostname, state->lhs_rank_world,
      redset_hostname, redset_rank,
      state->rhs_hostname, state->rhs_rank_world
    );
  }

  return rc;
}

#define REDSET_KEY_COPY_XOR_RANKS "RANKS"
#define REDSET_KEY_COPY_XOR_GROUP "GROUP"
#define REDSET_KEY_COPY_XOR_GROUP_RANK  "RANK"
#define REDSET_KEY_COPY_XOR_GROUP_RANKS "RANKS"

/* given a redundancy descriptor with all top level fields filled in
 * allocate and fill in structure for xor specific fields in state */
static int redset_create_xor(MPI_Comm parent_comm, redset* d)
{
  int rc = REDSET_SUCCESS;

  /* allocate a new structure to hold XOR state */
  redset_xor* state = (redset_xor*) REDSET_MALLOC(sizeof(redset_xor));

  /* attach structure to reddesc */
  d->state = (void*) state;

  /* allocate a new hash to store group mapping info */
  kvtree* header = kvtree_new();

  /* record the total number of ranks in parent comm */
  int parent_ranks;
  MPI_Comm_size(parent_comm, &parent_ranks);
  kvtree_set_kv_int(header, REDSET_KEY_COPY_XOR_RANKS, parent_ranks);

  /* create a new empty hash to track group info for this xor set */
  kvtree* hash = kvtree_new();
  kvtree_set(header, REDSET_KEY_COPY_XOR_GROUP, hash);

  /* record the total number of ranks in the xor communicator */
  int ranks_comm;
  MPI_Comm_size(d->comm, &ranks_comm);
  kvtree_set_kv_int(hash, REDSET_KEY_COPY_XOR_GROUP_RANKS, ranks_comm);

  /* record mapping of rank in xor group to corresponding parent rank */
  if (ranks_comm > 0) {
    /* allocate array to receive rank from each process */
    int* ranklist = (int*) REDSET_MALLOC(ranks_comm * sizeof(int));

    /* gather rank values */
    int parent_rank;
    MPI_Comm_rank(parent_comm, &parent_rank);
    MPI_Allgather(&parent_rank, 1, MPI_INT, ranklist, 1, MPI_INT, d->comm);

    /* map ranks in comm to ranks in comm */
    int i;
    for (i=0; i < ranks_comm; i++) {
      int rank = ranklist[i];
      kvtree_setf(hash, NULL, "%s %d %d", REDSET_KEY_COPY_XOR_GROUP_RANK, i, rank);
    }

    /* free the temporary array */
    redset_free(&ranklist);
  }

  /* record group mapping info in descriptor */
  state->group_map = header; 

  /* record group rank, world rank, and hostname of left and right partners */
  redset_set_partners(
    parent_comm, d->comm, 1,
    &state->lhs_rank, &state->lhs_rank_world, &state->lhs_hostname,
    &state->rhs_rank, &state->rhs_rank_world, &state->rhs_hostname
  );

  /* check that we got valid partners */
  if (state->lhs_hostname == NULL ||
      state->rhs_hostname == NULL ||
      strcmp(state->lhs_hostname, "") == 0 ||
      strcmp(state->rhs_hostname, "") == 0 ||
      strcmp(state->lhs_hostname, redset_hostname) == 0 ||
      strcmp(state->rhs_hostname, redset_hostname) == 0)
  {
    /* disable this descriptor */
    d->enabled = 0;
    redset_warn("Failed to find partner processes for redundancy descriptor, disabling, too few nodes? @ %s:%d",
      __FILE__, __LINE__
    );
    rc = REDSET_FAILURE;
  } else {
    redset_dbg(2, "LHS partner: %s (%d)  -->  My name: %s (%d)  -->  RHS partner: %s (%d)",
      state->lhs_hostname, state->lhs_rank_world,
      redset_hostname, redset_rank,
      state->rhs_hostname, state->rhs_rank_world
    );
  }

  return rc;
}

static int redset_delete_partner(redset* d)
{
  redset_partner* state = (redset_partner*) d->state;
  if (state != NULL) {
    /* free strings that we received */
    redset_free(&state->lhs_hostname);
    redset_free(&state->rhs_hostname);

    /* free the structure */
    redset_free(&d->state);
  }
  return REDSET_SUCCESS;
}

static int redset_delete_xor(redset* d)
{
  redset_xor* state = (redset_xor*) d->state;
  if (state != NULL) {
    /* free the hash mapping group ranks to world ranks */
    kvtree_delete(&state->group_map);

    /* free strings that we received */
    redset_free(&state->lhs_hostname);
    redset_free(&state->rhs_hostname);

    /* free the structure */
    redset_free(&d->state);
  }
  return REDSET_SUCCESS;
}

/* free any memory associated with the specified redundancy
 * descriptor */
int redset_delete(redset* d)
{
  /* free off copy type specific data */
  switch (d->type) {
  case REDSET_COPY_SINGLE:
    break;
  case REDSET_COPY_PARTNER:
    redset_delete_partner(d);
    break;
  case REDSET_COPY_XOR:
    redset_delete_xor(d);
    break;
  }

  /* free the communicator we created */
  if (d->comm != MPI_COMM_NULL) {
    MPI_Comm_free(&d->comm);
  }

  return REDSET_SUCCESS;
}

#define REDSET_KEY_CONFIG_ENABLED   "ENABLED"
#define REDSET_KEY_CONFIG_INTERVAL  "INTERVAL"
#define REDSET_KEY_CONFIG_OUTPUT    "OUTPUT"
#define REDSET_KEY_CONFIG_STORE     "STORE"
#define REDSET_KEY_CONFIG_DIRECTORY "DIR"
#define REDSET_KEY_CONFIG_TYPE      "TYPE"
#define REDSET_KEY_CONFIG_SET_SIZE  "SETSIZE"
#define REDSET_KEY_CONFIG_GROUP      "GROUP"
#define REDSET_KEY_CONFIG_GROUPS     "GROUPS"
#define REDSET_KEY_CONFIG_GROUP_ID   "GROUP"
#define REDSET_KEY_CONFIG_GROUP_SIZE "RANKS"
#define REDSET_KEY_CONFIG_GROUP_RANK "RANK"

/* convert the specified redundancy descritpor into a corresponding
 * kvtree */
int redset_to_kvtree(const redset* d, kvtree* hash)
{
  /* check that we got a valid pointer to a redundancy descriptor and
   * a hash */
  if (d == NULL || hash == NULL) {
    return REDSET_FAILURE;
  }

  /* clear the hash */
  kvtree_unset_all(hash);

  /* set the ENABLED key */
  kvtree_set_kv_int(hash, REDSET_KEY_CONFIG_ENABLED, d->enabled);

  /* set the TYPE key */
  switch (d->type) {
  case REDSET_COPY_SINGLE:
    kvtree_set_kv(hash, REDSET_KEY_CONFIG_TYPE, "SINGLE");
    break;
  case REDSET_COPY_PARTNER:
    kvtree_set_kv(hash, REDSET_KEY_CONFIG_TYPE, "PARTNER");
    break;
  case REDSET_COPY_XOR:
    kvtree_set_kv(hash, REDSET_KEY_CONFIG_TYPE, "XOR");
    break;
  }

  /* we don't set the LHS or RHS values because they are dependent on
   * runtime environment */

  /* we don't set the COMM because this is dependent on runtime
   * environment */

  /* set the GROUP_ID and GROUP_RANK keys, we use this info to rebuild
   * our communicator later */
  kvtree_set_kv_int(hash, REDSET_KEY_CONFIG_GROUPS,     d->groups);
  kvtree_set_kv_int(hash, REDSET_KEY_CONFIG_GROUP_ID,   d->group_id);
  kvtree_set_kv_int(hash, REDSET_KEY_CONFIG_GROUP_SIZE, d->ranks);
  kvtree_set_kv_int(hash, REDSET_KEY_CONFIG_GROUP_RANK, d->rank);

  return REDSET_SUCCESS;
}

/* given our rank within a set of ranks and minimum group size,
 * divide the set as evenly as possible and return the
 * group id corresponding to our rank */
static int redset_group_id(
  int rank,
  int ranks,
  int minsize,
  int* group_id)
{
  /* compute maximum number of full minsize groups we can fit within
   * ranks */
  int groups = ranks / minsize;

  /* compute number of ranks left over */
  int remainder_ranks = ranks - groups * minsize;

  /* determine base size for each group */
  int size = ranks;
  if (groups > 0) {
    /* evenly distribute remaining ranks over the groups that we have */
    int add_to_each_group = remainder_ranks / groups;
    size = minsize + add_to_each_group;
  }

  /* compute remaining ranks assuming we have groups of the new base
   * size */
  int remainder = ranks % size;

  /* for each remainder rank, we increase the lower groups by a size
   * of one, so that we end up with remainder groups of size+1 followed
   * by (groups - remainder) of size */

  /* cutoff is the first rank for which all groups are exactly size */
  int cutoff = remainder * (size + 1);

  if (rank < cutoff) {
    /* ranks below cutoff are grouped into sets of size+1 */
    *group_id = rank / (size + 1);
  } else {
    /* ranks at cutoff and higher are grouped into sets of size */
    *group_id = (rank - cutoff) / size + remainder;
  }

  return REDSET_SUCCESS;
}

/* given a parent communicator and a communicator representing our group
 * within the parent, split parent into other communicators consisting
 * of all procs with same rank within its group */
static int redset_split_across(
  MPI_Comm comm_parent,
  MPI_Comm comm_group,
  MPI_Comm* comm_across)
{
  /* TODO: this works well if each comm has about the same number of
   * procs, but we need something better to handle unbalanced groups */

  /* get rank of this process within parent communicator */
  int rank_parent;
  MPI_Comm_rank(comm_parent, &rank_parent);

  /* get rank of this process within group communicator */
  int rank_group;
  MPI_Comm_rank(comm_group, &rank_group);

  /* Split procs in parent into groups containing all procs with same
   * rank within group, order by rank in parent */
  MPI_Comm_split(comm_parent, rank_group, rank_parent, comm_across);
  
  return REDSET_SUCCESS;
}

/* convert copy type string to integer value */
static int redset_type_int_from_str(const char* value, int* outtype)
{
  int rc = REDSET_SUCCESS;

  int type;
  if (strcasecmp(value, "SINGLE") == 0) {
    type = REDSET_COPY_SINGLE;
  } else if (strcasecmp(value, "PARTNER") == 0) {
    type = REDSET_COPY_PARTNER;
  } else if (strcasecmp(value, "XOR") == 0) {
    type = REDSET_COPY_XOR;
  } else {
    if (redset_rank == 0) {
      redset_warn("Unknown copy type %s @ %s:%d",
        value, __FILE__, __LINE__
      );
    }
    rc = REDSET_FAILURE;
  }

  *outtype = type;
  return rc;
}

/* build a redundancy descriptor corresponding to the specified hash,
 * this function is collective */
int redset_create(
  int type,
  MPI_Comm comm,
  const char* group_name,
  redset* d)
{
  int SET_SIZE = 8;

  int rc = REDSET_SUCCESS;

  /* initialize the descriptor */
  redset_initialize(d);

  /* set xor set size */
  int set_size = SET_SIZE;

  /* assume it's enabled, we may turn this bit off later */
  d->enabled = 1;

  /* record type of descriptor */
  d->type = type;

  /* split procs from comm into sub communicators based on group_name */
  MPI_Comm newcomm;
  int groups, group_id;
  mpi_rank_str(comm, group_name, &groups, &group_id);
  MPI_Comm_split(comm, group_id, 0, &newcomm);

  /* build the communicator based on the copy type
   * and other parameters */
  MPI_Comm comm_across;
  int rank_across, ranks_across, split_id;
  switch (d->type) {
  case REDSET_COPY_SINGLE:
    /* not going to communicate with anyone, so just dup COMM_SELF */
    MPI_Comm_dup(MPI_COMM_SELF, &d->comm);
    break;
  case REDSET_COPY_PARTNER:
    /* dup the communicator across failure groups */
    redset_split_across(comm, newcomm, &d->comm);
    break;
  case REDSET_COPY_XOR:
    /* split the communicator across nodes based on xor set size
     * to create our xor communicator */
    /* split comm world across failure groups */
    redset_split_across(comm, newcomm, &comm_across);

    /* get our rank and the number of ranks in this communicator */
    MPI_Comm_rank(comm_across, &rank_across);
    MPI_Comm_size(comm_across, &ranks_across);

    /* identify which group we're in */
    redset_group_id(rank_across, ranks_across, set_size, &split_id);

    /* split communicator into sets */
    MPI_Comm_split(comm_across, split_id, redset_rank, &d->comm);

    /* free the temporary communicator */
    MPI_Comm_free(&comm_across);
    break;
  }

  /* find our position in the reddesc communicator */
  MPI_Comm_rank(d->comm, &d->rank);
  MPI_Comm_size(d->comm, &d->ranks);

  /* for our group id, use the global rank of the rank 0 task
   * in our reddesc comm */
  int comm_rank;
  MPI_Comm_rank(comm, &comm_rank);
  d->group_id = comm_rank;
  MPI_Bcast(&d->group_id, 1, MPI_INT, 0, d->comm);

  /* count the number of groups */
  int group_master = (d->rank == 0) ? 1 : 0;
  MPI_Allreduce(&group_master, &d->groups, 1, MPI_INT, MPI_SUM, comm);

  /* fill in state struct depending on copy type */
  switch (d->type) {
  case REDSET_COPY_SINGLE:
    break;
  case REDSET_COPY_PARTNER:
    redset_create_partner(comm, d);
    break;
  case REDSET_COPY_XOR:
    redset_create_xor(comm, d);
    break;
  }

  /* free communicator of all procs in the same group */
  MPI_Comm_free(&newcomm);

  return REDSET_SUCCESS;
}

/* convert the specified redundancy descritpor into a corresponding
 * kvtree */
int redset_store_to_kvtree(const redset* d, kvtree* hash)
{
  /* check that we got a valid pointer to a redundancy descriptor and
   * a hash */
  if (d == NULL || hash == NULL) {
    return REDSET_FAILURE;
  }

  /* clear the hash */
  kvtree_unset_all(hash);

  /* set the ENABLED key */
  kvtree_set_kv_int(hash, REDSET_KEY_CONFIG_ENABLED, d->enabled);

  /* set the TYPE key */
  switch (d->type) {
  case REDSET_COPY_SINGLE:
    kvtree_set_kv(hash, REDSET_KEY_CONFIG_TYPE, "SINGLE");
    break;
  case REDSET_COPY_PARTNER:
    kvtree_set_kv(hash, REDSET_KEY_CONFIG_TYPE, "PARTNER");
    break;
  case REDSET_COPY_XOR:
    kvtree_set_kv(hash, REDSET_KEY_CONFIG_TYPE, "XOR");
    break;
  }

  /* we don't set the LHS or RHS values because they are dependent on
   * runtime environment */

  /* we don't set the COMM because this is dependent on runtime
   * environment */

  /* set the GROUP_ID and GROUP_RANK keys, we use this info to rebuild
   * our communicator later */
  kvtree_set_kv_int(hash, REDSET_KEY_CONFIG_GROUPS,     d->groups);
  kvtree_set_kv_int(hash, REDSET_KEY_CONFIG_GROUP_ID,   d->group_id);
  kvtree_set_kv_int(hash, REDSET_KEY_CONFIG_GROUP_SIZE, d->ranks);
  kvtree_set_kv_int(hash, REDSET_KEY_CONFIG_GROUP_RANK, d->rank);

  return REDSET_SUCCESS;
}

/* build a redundancy descriptor corresponding to the specified kvtree,
 * this function is collective, it is the opposite of store_to_kvtree */
int redset_restore_from_kvtree(
  MPI_Comm comm,
  const kvtree* hash,
  redset* d)
{
  int rc = REDSET_SUCCESS;

  /* check that we got a valid redundancy descriptor */
  if (d == NULL) {
    printf("No redundancy descriptor to fill from hash @ %s:%d",
      __FILE__, __LINE__
    );
    rc = REDSET_FAILURE;
  }

  /* check that we got a valid pointer to a hash */
  if (hash == NULL) {
    printf("No hash specified to build redundancy descriptor from @ %s:%d",
      __FILE__, __LINE__
    );
    rc = REDSET_FAILURE;
  }

  /* check that everyone made it this far */
  if (! redset_alltrue(rc == REDSET_SUCCESS, comm)) {
    if (d != NULL) {
      d->enabled = 0;
    }
    return REDSET_FAILURE;
  }

  /* initialize the descriptor */
  redset_initialize(d);

  /* enable / disable the descriptor */
  d->enabled = 1;
  kvtree_util_get_int(hash, REDSET_KEY_CONFIG_ENABLED, &(d->enabled));

  /* read the redundancy scheme type from the hash,
   * and build our redundancy communicator */
  char* type;
  if (kvtree_util_get_str(hash, REDSET_KEY_CONFIG_TYPE, &type) == REDSET_SUCCESS) {
    if (redset_type_int_from_str(type, &d->type) != REDSET_SUCCESS) {
      d->enabled = 0;
      if (redset_rank == 0) {
        redset_abort(-1, "Unknown copy type %s in redundancy descriptor hash @ %s:%d",
          type, __FILE__, __LINE__
        );
      }
    }
  } else {
    redset_abort(-1, "Unknown copy type in redundancy descriptor hash @ %s:%d",
      __FILE__, __LINE__
    );
  }

  /* build the group communicator */
  if (kvtree_util_get_int(hash, REDSET_KEY_CONFIG_GROUPS, &d->groups) != REDSET_SUCCESS) {
    redset_abort(-1, "Failed to read group count in redundancy descriptor hash @ %s:%d",
      __FILE__, __LINE__
    );
  }
  if (kvtree_util_get_int(hash, REDSET_KEY_CONFIG_GROUP_ID, &d->group_id) != REDSET_SUCCESS) {
    redset_abort(-1, "Failed to read group id in redundancy descriptor hash @ %s:%d",
      __FILE__, __LINE__
    );
  }
  if (kvtree_util_get_int(hash, REDSET_KEY_CONFIG_GROUP_RANK, &d->rank) != REDSET_SUCCESS) {
    redset_abort(-1, "Failed to read group rank in redundancy descriptor hash @ %s:%d",
      __FILE__, __LINE__
    );
  }
  if (kvtree_util_get_int(hash, REDSET_KEY_CONFIG_GROUP_SIZE, &d->ranks) != REDSET_SUCCESS) {
    redset_abort(-1, "Failed to read group size in redundancy descriptor hash @ %s:%d",
      __FILE__, __LINE__
    );
  }
  MPI_Comm_split(comm, d->group_id, d->rank, &d->comm);

  /* fill in state struct depending on copy type */
  switch (d->type) {
  case REDSET_COPY_SINGLE:
    break;
  case REDSET_COPY_PARTNER:
    redset_create_partner(comm, d);
    break;
  case REDSET_COPY_XOR:
    redset_create_xor(comm, d);
    break;
  }

  /* if anyone has disabled this, everyone needs to */
  if (! redset_alltrue(d->enabled, comm)) {
    d->enabled = 0;
  }

  return REDSET_SUCCESS;
}

/* checks whether specifed file exists, is readable, and is complete */
int redset_bool_have_file(const char* file)
{
  /* if no filename is given return false */
  if (file == NULL || strcmp(file,"") == 0) {
    redset_dbg(2, "File name is null or the empty string @ %s:%d",
      __FILE__, __LINE__
    );
    return 0;
  }

  /* check that we can read the file */
  if (redset_file_is_readable(file) != REDSET_SUCCESS) {
    redset_dbg(2, "Do not have read access to file: %s @ %s:%d",
      file, __FILE__, __LINE__
    );
    return 0;
  }

#if 0
  /* check that the file size matches */
  unsigned long size = scr_file_size(file);
  unsigned long meta_size = 0;
  if (scr_meta_get_filesize(meta, &meta_size) != SCR_SUCCESS) {
    scr_dbg(2, "Failed to read filesize field in meta data: %s @ %s:%d",
      file, __FILE__, __LINE__
    );
    scr_meta_delete(&meta);
    return 0;
  }
  if (size != meta_size) {
    scr_dbg(2, "Filesize is incorrect, currently %lu, expected %lu for %s @ %s:%d",
      size, meta_size, file, __FILE__, __LINE__
    );
    scr_meta_delete(&meta);
    return 0;
  }

  /* TODO: check that crc32 match if set (this would be expensive) */
#endif

  /* if we made it here, assume the file is good */
  return 1;
}

/* check whether we have all files for a given rank of a given dataset */
int redset_bool_have_files(
  int numfiles,
  const char** files,
  const char* name,
  const redset* d)
{
  /* check the integrity of each of the files */
  int missing_a_file = 0;
  int i;
  for (i = 0; i < numfiles; i++) {
    const char* file = files[i];
    if (! redset_bool_have_file(file)) {
      missing_a_file = 1;
    }
  }
  if (missing_a_file) {
    return 0;
  }

  /* if we make it here, we have all of our files */
  return 1;
}

static int redset_build_filename(
  const char* name,
  char* file,
  size_t len)
{
  snprintf(file, len, "%s.redset", name);
  return REDSET_SUCCESS;
}

/* encode redundancy descriptor and write to a file */
static int redset_encode_reddesc(
  const char* name,
  const redset* d)
{
  /* get name of this process */
  int rank_world;
  MPI_Comm comm_world = redset_comm;
  MPI_Comm_rank(comm_world, &rank_world);

  /* allocate a structure to record meta data about our files and redundancy descriptor */
  kvtree* current_hash = kvtree_new();

  /* store our redundancy descriptor in hash */
  kvtree* desc_hash = kvtree_new();
  redset_store_to_kvtree(d, desc_hash);
  kvtree_set(current_hash, "DESC", desc_hash);

  /* copy meta data to hash */
  kvtree* meta_hash = kvtree_new();
  kvtree_setf(meta_hash, current_hash, "%d", rank_world);

  /* apply redundancy to hash data according to selected scheme */
  int rc = REDSET_FAILURE;
  switch (d->type) {
  case REDSET_COPY_SINGLE:
    rc = redset_encode_reddesc_single(meta_hash, name, d);
    break;
  case REDSET_COPY_PARTNER:
    rc = redset_encode_reddesc_partner(meta_hash, name, d);
    break;
  case REDSET_COPY_XOR:
    rc = redset_encode_reddesc_xor(meta_hash, name, d);
    break;
  }

  /* write meta data to file */
  char filename[REDSET_MAX_FILENAME];
  redset_build_filename(name, filename, sizeof(filename));
  kvtree_write_file(filename, meta_hash);

  /* delete the hash */
  kvtree_delete(&meta_hash);

  return rc;
}

/* delete files added in encode_reddesc */
static int redset_unencode_reddesc(
  const char* name,
  const redset* d)
{
  /* delete meta data file */
  char filename[REDSET_MAX_FILENAME];
  redset_build_filename(name, filename, sizeof(filename));
  int rc = redset_file_unlink(filename);
  return rc;
}

/* attempt to recover and rebuild redundancy descriptor information */
static int redset_recover_reddesc(
  const char* name,
  redset* d)
{
  /* get name of this process */
  int rank_world;
  MPI_Comm comm_world = redset_comm;
  MPI_Comm_rank(comm_world, &rank_world);

  /* read meta data from file */
  kvtree* hash = kvtree_new();
  char filename[REDSET_MAX_FILENAME];
  redset_build_filename(name, filename, sizeof(filename));
  kvtree_read_file(filename, hash);

  /* create a hash to exchange redundancy descriptors */
  kvtree* send_hash = kvtree_new();

  /* determine which processes we have redundancy descriptors for,
   * and include them in the send hash */
  kvtree_elem* elem;
  for (elem = kvtree_elem_first(hash);
       elem != NULL;
       elem = kvtree_elem_next(elem))
  {
    /* get rank of this element */
    int rank = kvtree_elem_key_int(elem);

    kvtree* rank_hash = kvtree_elem_hash(elem);

    /* see if we have a descriptor for this rank */
    kvtree* desc_hash = kvtree_get(rank_hash, "DESC");
    if (desc_hash != NULL) {
      kvtree* desc_copy = kvtree_new();
      kvtree_merge(desc_copy, desc_hash);
      kvtree_setf(send_hash, desc_copy, "%d", rank);
//      kvtree_exchange_sendq(send_hash, rank, desc_hash);
    }
  }

  kvtree_delete(&hash);

  /* create an empty hash to receive any incoming descriptors */
  /* exchange descriptors with other ranks */
  kvtree* recv_hash = kvtree_new();
  kvtree_exchange(send_hash, recv_hash, comm_world);

  /* check that everyone can get their descriptor */
  int num_desc = kvtree_size(recv_hash);
  if (! redset_alltrue((num_desc > 0), comm_world)) {
    kvtree_delete(&recv_hash);
    kvtree_delete(&send_hash);
    redset_dbg(2, "Cannot find process that has my redundancy descriptor @ %s:%d",
      __FILE__, __LINE__
    );
    return REDSET_FAILURE;
  }

  /* just go with the first redundancy descriptor in our list,
   * they should all be the same */
  kvtree_elem* desc_elem = kvtree_elem_first(recv_hash);
  kvtree* desc_hash = kvtree_elem_hash(desc_elem);

  /* rebuild the redundancy descriptor for this dataset */
  redset_restore_from_kvtree(comm_world, desc_hash, d);

  /* can free our hashes now that we've rebuilt our redundancy descriptor */
  kvtree_delete(&recv_hash);
  kvtree_delete(&send_hash);

  /* reapply encoding to redundancy descriptor */
  int rc = redset_encode_reddesc(name, d);

  return rc;
}

/* apply redundancy scheme to list of files and given redundancy
 * descriptor, scheme may store redundany data
 * and associate it with given name */
int redset_apply(
  int numfiles,
  const char** files,
  const char* name,
  const redset* d)
{
  int i;
  MPI_Comm comm_world = redset_comm;

  int rank_world, nranks_world;
  MPI_Comm_rank(comm_world, &rank_world);
  MPI_Comm_size(comm_world, &nranks_world);

  /* start timer */
  time_t timestamp_start;
  double time_start;
  if (rank_world == 0) {
    time_start = MPI_Wtime();
  }

  /* initialize to 0 */
  double bytes = 0.0;

  /* step through each of my files for the specified dataset
   * to scan for any incomplete files */
  int valid = 1;
  double my_bytes = 0.0;
  for (i = 0; i < numfiles; i++) {
    /* get file name of this file */
    const char* file_name = files[i];

    /* TODO: check that file exists */

    /* get file size of this file */
    unsigned long file_size = redset_file_size(file_name);

    /* add up the number of bytes on our way through */
    my_bytes += (double) file_size;
  }

  /* TODO: determine whether everyone's files are good */

  /* apply encoding to redundancy descriptor */
  int rc = redset_encode_reddesc(name, d);

  /* determine whether all processes saved their redundancy info */
  if (! redset_alltrue((rc == REDSET_SUCCESS), comm_world)) {
    /* at least one process failed to rebuild its redundancy information */
    return REDSET_FAILURE;
  }

  /* apply the redundancy scheme */
  switch (d->type) {
  case REDSET_COPY_SINGLE:
    rc = redset_apply_single(numfiles, files, name, d);
    break;
  case REDSET_COPY_PARTNER:
    rc = redset_apply_partner(numfiles, files, name, d);
    break;
  case REDSET_COPY_XOR:
    rc = redset_apply_xor(numfiles, files, name, d);
    break;
  }

  /* determine whether everyone succeeded in their copy */
  if (! redset_alltrue((rc == REDSET_SUCCESS), comm_world)) {
    /* at least one process failed to apply redundancy scheme */
    return rc = REDSET_FAILURE;
  }

  /* add up total number of bytes */
  MPI_Allreduce(&my_bytes, &bytes, 1, MPI_DOUBLE, MPI_SUM, comm_world);

  /* stop timer and report performance info */
  if (rank_world == 0) {
    double time_end = MPI_Wtime();
    double time_diff = time_end - time_start;
    double bw = 0.0;
    if (time_diff > 0) {
      bytes / (1024.0 * 1024.0 * time_diff);
    }
    //printf("redset_apply: %f secs, %e bytes, %f MB/s, %f MB/s per proc\n",
    //        time_diff, bytes, bw, bw/(double)nranks_world
    //);

    // TODO: log or report cost somewhere
  }

  return rc;
}

/* rebuilds files for specified dataset, return redundancy descriptor,
 * returns REDSET_SUCCESS if successful, REDSET_FAILURE otherwise,
 * the same return code is returned to all processes */
int redset_recover(
  const char* name,
  redset* d)
{
  MPI_Comm comm_world = redset_comm;

  /* reapply encoding to redundancy descriptor */
  int rc = redset_recover_reddesc(name, d);

  /* determine whether everyone succeeded */
  if (! redset_alltrue((rc == REDSET_SUCCESS), comm_world)) {
    /* at least one process failed to rebuild its redundancy information */
    if (rc == REDSET_SUCCESS) {
      /* free our descriptor if we were successful */
      redset_delete(d);
    }
    return REDSET_FAILURE;
  }

  /* now recover data using redundancy scheme, if necessary */
  switch (d->type) {
  case REDSET_COPY_SINGLE:
    rc = redset_recover_single(name, d);
    break;
  case REDSET_COPY_PARTNER:
    rc = redset_recover_partner(name, d);
    break;
  case REDSET_COPY_XOR:
    rc = redset_recover_xor(name, d);
    break;
  }

  /* determine whether everyone succeeded */
  if (! redset_alltrue((rc == REDSET_SUCCESS), comm_world)) {
    /* at least one process failed to rebuild its data */
    redset_delete(d);
    return REDSET_FAILURE;
  }

  return REDSET_SUCCESS;
}

/* deletes redundancy data that was added in redset_apply,
 * which is useful when cleaning up */
int redset_unapply(
  const char* name,
  redset* d)
{
  int rc;
  MPI_Comm comm_world = redset_comm;

  /* now recover data using redundancy scheme, if necessary */
  switch (d->type) {
  case REDSET_COPY_SINGLE:
    rc = redset_unapply_single(name, d);
    break;
  case REDSET_COPY_PARTNER:
    rc = redset_unapply_partner(name, d);
    break;
  case REDSET_COPY_XOR:
    rc = redset_unapply_xor(name, d);
    break;
  }

  /* determine whether everyone succeeded */
  if (! redset_alltrue((rc == REDSET_SUCCESS), comm_world)) {
    /* at least one process failed to rebuild its data */
    return REDSET_FAILURE;
  }

  /* reapply encoding to redundancy descriptor */
  rc = redset_unencode_reddesc(name, d);

  /* determine whether everyone succeeded */
  if (! redset_alltrue((rc == REDSET_SUCCESS), comm_world)) {
    /* at least one process failed to rebuild its redundancy information */
    return REDSET_FAILURE;
  }

  return REDSET_SUCCESS;
}

static int redset_from_dir(
  const char* name, 
  redset* d)
{
  /* get name of this process */
  int rank;
  MPI_Comm comm_world = redset_comm;
  MPI_Comm_rank(comm_world, &rank);

  char rankstr[100];
  snprintf(rankstr, sizeof(rankstr), "%d", rank);

  /* read meta data from file */
  kvtree* hash = kvtree_new();
  char filename[REDSET_MAX_FILENAME];
  redset_build_filename(name, filename, sizeof(filename));
  kvtree_read_file(filename, hash);

  /* get pointer to descriptor hash for this rank */
  kvtree* desc_hash = kvtree_get_kv(hash, rankstr, "DESC");

  /* rebuild the redundancy descriptor for this dataset */
  int rc = redset_restore_from_kvtree(comm_world, desc_hash, d);

  /* free the hash */
  kvtree_delete(&hash);

  return rc;
}

/* returns a list of files added by redundancy descriptor */
redset_filelist* redset_filelist_get(
  const char* name,
  redset* d)
{
  /* build redundancy descriptor from name */
  //redset d;
  //redset_from_dir(name, &d);

  /* create a temporary file list to record files from redundancy scheme */
  redset_filelist* tmp = NULL;

  /* get files added by redundancy method */
  switch (d->type) {
  case REDSET_COPY_SINGLE:
    tmp = redset_filelist_get_single(name, d);
    break;
  case REDSET_COPY_PARTNER:
    tmp = redset_filelist_get_partner(name, d);
    break;
  case REDSET_COPY_XOR:
    tmp = redset_filelist_get_xor(name, d);
    break;
  }

  /* we have a top level redundancy file in addition to anything added by scheme,
   * allocate space for the full file list */
  int count = tmp->count + 1;
  const char** files = (const char**) REDSET_MALLOC(count * sizeof(char*));

  /* record name of top level file */
  char filename[REDSET_MAX_FILENAME];
  redset_build_filename(name, filename, sizeof(filename));
  files[0] = strdup(filename);

  /* record each redundancy file */
  int i;
  for (i = 0; i < tmp->count; i++) {
    files[i+1] = strdup(tmp->files[i]);
  }

  /* free the temporary list */
  redset_filelist_release(&tmp);

  /* free the descriptor */
  //redset_free(&d);

  /* update fields of return list */
  redset_filelist* list = (redset_filelist*) REDSET_MALLOC(sizeof(redset_filelist));
  list->count = count;
  list->files = files;

  return list;
}

int redset_filelist_release(redset_filelist** plist)
{
  if (plist == NULL) {
    return REDSET_SUCCESS;
  }

  redset_filelist* list = *plist;

  /* check that we got a list */
  if (list == NULL) {
    return REDSET_SUCCESS;
  }

  /* free each file name string */
  int i;
  for (i = 0; i < list->count; i++) {
    redset_free(&list->files[i]);
  }

  /* free the list of files itself */
  redset_free(&list->files);

  /* free the object */
  redset_free(plist);

  return REDSET_SUCCESS;
}

/* returns the number of files in the list */
int redset_filelist_count(redset_filelist* list)
{
  /* check that we got a list */
  if (list == NULL) {
    return 0;
  }

  return list->count;
}

/* returns the name of the file by the given index,
 * index should be between 0 and count-1,
 * returns NULL if index is invalid */
const char* redset_filelist_file(redset_filelist* list, int index)
{
  /* check that we got a list */
  if (list == NULL) {
    return NULL;
  }

  /* check that index is in range */
  if (index < 0 || index >= list->count) {
    return NULL;
  }

  return list->files[index];
}
