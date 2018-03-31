#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "kvtree.h"
#include "kvtree_util.h"
#include "mpi.h"
#include "redset.h"
#include "redset_util.h"

static int er_scheme_counter = 0;
static int er_set_counter = 0;

static kvtree* er_schemes = NULL;
static kvtree* er_sets = NULL;

int ER_Init(const char* conf_file)
{
  int rc = ER_SUCCESS;

  /* allocate maps to track scheme and set data */
  er_schemes = kvtree_new();
  er_sets    = kvtree_new();

  /* initialize the redundancy library */
  if (redset_init() != REDSET_SUCCESS) {
    rc = ER_FAILURE;
  }

  return rc;
}

int ER_Finalize()
{
  int rc = ER_SUCCESS;

  /* TODO: free descriptors for any outstanding schemes,
   * probably need to do this in same order on all procs,
   * for now, force user to clean up */
  kvtree* schemes = kvtree_get(er_schemes, "SCHEMES");
  if (kvtree_size(schemes) > 0) {
    redset_err("ER_Finalize called before schemes freed");
    rc = ER_FAILURE;
  }

  /* free maps */
  kvtree_delete(&er_schemes);
  kvtree_delete(&er_sets);

  /* shut down redundancy library */
  if (redset_finalize() != REDSET_SUCCESS) {
    rc = ER_FAILURE;
  }

  return rc;
}

int ER_Create_Scheme(
  MPI_Comm comm,
  const char* failure_domain,
  int encoding_blocks,
  int erasure_blocks)
{
  int rc = ER_SUCCESS;

  /* check that we can support the scheme the caller is asking for */
  int encoding_type = REDSET_COPY_NULL;
  if (encoding_blocks < 1) {
    /* no data to be encoded, don't know what to do */
    return ER_FAILURE;
  }
  if (erasure_blocks == 0) {
    encoding_type = REDSET_COPY_SINGLE;
  } else if (encoding_blocks == erasure_blocks) {
    encoding_type = REDSET_COPY_PARTNER;
  } else if (erasure_blocks == 1) {
    encoding_type = REDSET_COPY_XOR;
  } else {
    /* some form of Reed-Solomon that we don't support yet */
    return ER_FAILURE;
  }

  /* allocate a new redundancy descriptor */
  redset* d = REDSET_MALLOC(sizeof(redset));

  /* create the scheme */
  if (redset_create(encoding_type, comm, failure_domain, d) == REDSET_SUCCESS) {
    /* bump our internal counter */
    er_scheme_counter++;

    /* create a new record in our scheme map */
    kvtree* scheme = kvtree_set_kv_int(er_schemes, "SCHEMES", er_scheme_counter);

    /* record pointer to reddesc in our map */
    kvtree_util_set_ptr(scheme, "PTR", (void*)d);

    rc = ER_SUCCESS;
  } else {
    rc = ER_FAILURE;
  }

  /* return scheme id to caller, -1 if error */
  int ret = er_scheme_counter;
  if (rc != ER_SUCCESS) {
    ret = -1;
  }
  return ret;
}

int ER_Free_Scheme(int scheme_id)
{
  int rc = ER_SUCCESS;

  /* look up entry for this scheme id */
  kvtree* scheme = kvtree_get_kv_int(er_schemes, "SCHEMES", scheme_id);
  if (scheme) {
    /* get pointer to reddesc */
    redset* d;
    if (kvtree_util_get_ptr(scheme, "PTR", (void**)&d) == KVTREE_SUCCESS) {
      /* free reddesc */
      if (redset_delete(d) != REDSET_SUCCESS) {
        /* failed to free redundancy descriptor */
        rc = ER_FAILURE;
      }
    } else {
      /* failed to find pointer to reddesc */
      rc = ER_FAILURE;
    }

    /* drop scheme entry from our map */
    kvtree_unset_kv_int(er_schemes, "SCHEMES", scheme_id);
  } else {
    /* failed to find scheme id in map */
    rc = ER_FAILURE;
  }

  return rc;
}

/* create a named set, and specify whether it should be encoded or recovered */
int ER_Create(const char* name, int direction, int scheme_id)
{
  int rc = ER_SUCCESS;

  /* check that we got a name */
  if (name == NULL || strcmp(name, "") == 0) {
    return ER_FAILURE;
  }

  /* check that we got a valid value for direction */
  if (direction != ER_DIRECTION_ENCODE &&
      direction != ER_DIRECTION_REBUILD)
  {
    return ER_FAILURE;
  }

  /* bump set counter */
  er_set_counter++;

  /* add an entry for this set */
  kvtree* set = kvtree_set_kv_int(er_sets, "SETS", er_set_counter);

  /* record operation direction */
  kvtree_util_set_str(set, "NAME", name);

  /* record operation direction */
  kvtree_util_set_int(set, "DIRECTION", direction);

  /* when encoding, we need to remember the scheme,
   * it's implied by name on rebuild */
  if (direction == ER_DIRECTION_ENCODE) {
    /* record scheme id */
    kvtree_util_set_int(set, "SCHEME", scheme_id);

    /* look up entry for this scheme id */
    kvtree* scheme = kvtree_get_kv_int(er_schemes, "SCHEMES", scheme_id);
    if (! scheme) {
      /* failed to find scheme id in map */
      rc = ER_FAILURE;
    }
  }

  /* return set id to caller, -1 if error */
  int ret = er_set_counter;
  if (rc != ER_SUCCESS) {
    ret = -1;
  }
  return ret;
}

/* adds file to specified set id */
int ER_Add(int set_id, const char* file)
{
  int rc = ER_SUCCESS;

  /* check that we got a file name */
  if (file == NULL || strcmp(file, "") == 0) {
    return ER_FAILURE;
  }

  /* lookup set id */
  kvtree* set = kvtree_get_kv_int(er_sets, "SETS", set_id);
  if (set) {
    /* add file to set */
    kvtree_set_kv(set, "FILE", file);

    /* TODO: capture current working dir? */
  } else {
    /* failed to find set id */
    rc = ER_FAILURE;
  }

  return rc;
}

/* initiate encode/rebuild operation on specified set id */
int ER_Dispatch(int set_id)
{
  int rc = ER_SUCCESS;

  /* lookup set id */
  kvtree* set = kvtree_get_kv_int(er_sets, "SETS", set_id);
  if (set) {
    /* get name of set */
    char* name = NULL;
    if (kvtree_util_get_str(set, "NAME", &name) != KVTREE_SUCCESS) {
      rc = ER_FAILURE;
    }

    /* TODO: allow caller to specify directory */
    char dir[1024];
    snprintf(dir, sizeof(dir), ".er.%s", name);

    /* get direction of operation (encode / rebuild) */
    int direction = 0;
    if (kvtree_util_get_int(set, "DIRECTION", &direction) != KVTREE_SUCCESS) {
      rc = ER_FAILURE;
    }

    if (direction == ER_DIRECTION_ENCODE) {
      /* determine number of files */
      kvtree* files_hash = kvtree_get(set, "FILE");
      int num_files = kvtree_size(files_hash);

      /* allocate space for file names */
      const char** filenames = (const char**) REDSET_MALLOC(num_files * sizeof(char*));

      /* copy pointers to filenames */
      int i = 0;
      kvtree_elem* elem;
      for (elem = kvtree_elem_first(files_hash);
           elem != NULL;
           elem = kvtree_elem_next(elem))
      {
        const char* file = kvtree_elem_key(elem);
        filenames[i] = file;
        i++;
      }

      /* get scheme id */
      int scheme_id = 0;
      if (kvtree_util_get_int(set, "SCHEME", &scheme_id) != KVTREE_SUCCESS) {
        rc = ER_FAILURE;
      }

      /* get scheme */
      kvtree* scheme = kvtree_get_kv_int(er_schemes, "SCHEMES", scheme_id);
      if (scheme) {
        /* get redundancy descriptor from scheme */
        redset* d = NULL;
        if (kvtree_util_get_ptr(scheme, "PTR", (void**)&d) != KVTREE_SUCCESS) {
          rc = ER_FAILURE;
        }

        if (rc == ER_SUCCESS) {
          /* make directory */
          mode_t mode_dir = redset_getmode(1, 1, 1);
          redset_mkdir(dir, mode_dir);

          /* TODO: read process name from scheme? */
          char proc_name[100];
          int rank;
          MPI_Comm_rank(MPI_COMM_WORLD, &rank);
          snprintf(proc_name, sizeof(proc_name), "%s/%d", dir, rank);

          /* apply redundancy */
          if (redset_apply(num_files, filenames, proc_name, *d) != REDSET_SUCCESS) {
            /* failed to apply redundancy descriptor */
            rc = ER_FAILURE;
          }
        }
      } else {
        /* failed to find scheme id for this set */
        rc = ER_FAILURE;
      }

      /* free list of file names */
      redset_free(&filenames);
    } else {
      /* TODO: read process name from scheme? */
      char proc_name[100];
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      snprintf(proc_name, sizeof(proc_name), "%s/%d", dir, rank);

      /* rebuild files */
      redset d;
      if (redset_recover(MPI_COMM_WORLD, proc_name, &d) != REDSET_SUCCESS) {
        /* rebuild failed */
        rc = ER_FAILURE;
      }
    }
  } else {
    /* failed to find set id */
    rc = ER_FAILURE;
  }

  return rc;
}

/* tests whether ongoing dispatch operation to finish,
 * returns 1 if done, 0 otherwise */
int ER_Test(int set_id)
{
  return 1;
}

/* wait for ongoing dispatch operation to finish */
int ER_Wait(int set_id)
{
  return ER_SUCCESS;
}

/* free internal resources associated with set id */
int ER_Free(int set_id)
{
  kvtree_unset_kv_int(er_sets, "SETS", set_id);
  return ER_SUCCESS;
}
