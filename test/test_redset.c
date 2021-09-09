#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include <string.h>
#include <unistd.h>

#include <limits.h>
#include <unistd.h>

/* compute crc32 */
#include <zlib.h>

#include "mpi.h"

#include "redset.h"
#include "redset_io.h"

/* print error message to stdout */
#define ABORT(...) test_abort(__FILE__, __LINE__, 1, __VA_ARGS__)
void test_abort(const char* file, int line, int rc, const char *fmt, ...)
{
  va_list argp;
  fprintf(stdout, "ERROR: ");
  va_start(argp, fmt);
  vfprintf(stdout, fmt, argp);
  va_end(argp);
  fprintf(stdout, " @ %s:%d\n", file, line);
  fflush(stdout);

  MPI_Abort(MPI_COMM_WORLD, rc);
}

/* compute MPI_LAND of value across all procs,
 * return 1 if value is 1 on all procs, 0 otherwise */
int alltrue(int value, MPI_Comm comm)
{
  int all_value;
  MPI_Allreduce(&value, &all_value, 1, MPI_INT, MPI_LAND, comm);
  return all_value;
}

int create_files(int count, const char** filelist)
{
  int rc = 0;

  int rank, ranks;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &ranks);

  srand(rank);
  size_t bufsize = 1024*1024;
  //size_t bufsize = 1;
  char* buf = (char*) malloc(bufsize);
  int j;
  for (j = 0; j < bufsize; j++) {
    buf[j] = (char) rand();
    //buf[j] = 'A' + (char)rank;
  }

  int i;
  for (i = 0; i < count; i++) {
    //char buf[256];
    //sprintf(buf, "data from rank %d in file %d\n", rank, i);

    const char* name = filelist[i];
    int fd = open(name, O_WRONLY | O_TRUNC | O_CREAT, S_IRUSR | S_IWUSR);
    if (fd != -1) {
      int k;
      for (k = 0; k < ranks + rank; k++) {
        //for (j = 0; j < bufsize; j++) {
        //  buf[j] = (char) rand();
        //}
        //write(fd, buf, strlen(buf));
        write(fd, buf, bufsize);
      }
      close(fd);
    } else {
      ABORT("opening file %s: errno=%d %s", name, errno, strerror(errno));
    }
  }

  free(buf);
  return 0;
}

void init_crcs(int count, const char** filelist, uLong* crcvals)
{
  int i;
  for (i = 0; i < count; i++) {
    const char* file = filelist[i];
    redset_crc32(file, &crcvals[i]);
  }
}

int check_crcs(int count, const char** filelist, uLong* crcvals)
{
  int i;
  for (i = 0; i < count; i++) {
    uLong c;
    const char* file = filelist[i];
    redset_crc32(file, &c);
    if (c != crcvals[i]) {
      ABORT("Failed CRC32 check for %s", file);
    }
  }
  return 0;
}

/* delete the list of files */
int delete_files(int count, const char** filelist)
{
  int i;
  for (i = 0; i < count; i++) {
    const char* name = filelist[i];
    unlink(name);
  }
  return 0;
}

/* check that we have each redundancy file */
int check_for_redundancy_files(int mode, const char* path, redset d)
{
  int rc = 0;

  /* get list of redundancy files */
  redset_filelist list = redset_filelist_get(path, d);
  if (list == NULL) {
    ABORT("Failed to get list of redundancy files");
  }

  /* check that each file in the list exists */
  int i;
  int count = redset_filelist_count(list);
  for (i = 0; i < count; i++) {
    const char* name = redset_filelist_file(list, i);
    if (access(name, F_OK) != 0) {
      ABORT("Missing redundancy file %s", name);
    }
  }

  /* check that we got the expected number of files in the list */
  if (count != 2) {
    /* all methods generate two files */
    ABORT("Unexpected number of redundancy files");
  }

  /* free the list */
  redset_filelist_release(&list);
  if (list != NULL) {
    ABORT("Failed to free list of redundancy files");
  }

  return rc;
}

int delete_redundancy_files(int mode, const char* path, redset d)
{
  int rc = 0;

  /* get list of redundancy files */
  redset_filelist list = redset_filelist_get(path, d);
  if (list == NULL) {
    ABORT("Failed to get list of redundancy files");
  }

  /* check that each file in the list exists */
  int i;
  int count = redset_filelist_count(list);
  for (i = 0; i < count; i++) {
    const char* name = redset_filelist_file(list, i);
    if (unlink(name) != 0) {
      ABORT("Deleting redundancy file %s", name);
    }
  }

  /* check that we got the expected number of files in the list */
  if (count != 2) {
    /* all methods generate two files */
    ABORT("Unexpected number of redundancy files");
  }

  /* free the list */
  redset_filelist_release(&list);
  if (list != NULL) {
    ABORT("Failed to free list of redundancy files");
  }

  return rc;
}

/* apply redundancy descriptor and check that redundancy files are created */
int test_apply(int mode, int k, int filecount, const char** filelist, uLong* crcvals, const char* path, redset d, MPI_Comm comm)
{
  int rc = 0;

  int rank;
  MPI_Comm_rank(comm, &rank);

  double start = MPI_Wtime();

  int redset_rc = redset_apply(filecount, filelist, path, d);
  if (redset_rc != REDSET_SUCCESS) {
    ABORT("apply failed");
  }

  double end = MPI_Wtime();
  if (rank == 0) {
    printf("mode = %d: Apply time (protect %d ranks): %f\n", mode, k, (end - start));
  }

  check_for_redundancy_files(mode, path, d);

  check_crcs(filecount, filelist, crcvals);

  return rc;
}

/* unapply redundancy descriptor and check that redundancy files are gone */
int test_unapply(int mode, const char* path, redset d, MPI_Comm comm)
{
  int rc = 0;

  int redset_rc = redset_unapply(path, d);
  if (! alltrue(redset_rc == REDSET_SUCCESS, comm)) {
    ABORT("unapply failed");
  }

  /* get list of redundancy files */
  redset_filelist list = redset_filelist_get(path, d);
  if (list == NULL) {
    ABORT("Failed to get list of redundancy files");
  }

  /* check that each file in the list is gone */
  int i;
  int count = redset_filelist_count(list);
  for (i = 0; i < count; i++) {
    const char* name = redset_filelist_file(list, i);
    if (access(name, F_OK) == 0) {
      ABORT("Found redundancy file %s", name);
    }
  }

  /* check that we got the expected number of files in the list */
  if (count != 2) {
    /* all methods generate two files */
    ABORT("Unexpected number of redundancy files");
  }

  /* free the list */
  redset_filelist_release(&list);
  if (list != NULL) {
    ABORT("Failed to free list of redundancy files");
  }

  return rc;
}

/* test that recover succeeds when no files are lost,
 * nothing to do, but should still return success */
int test_recover_no_loss(int mode, const char* path, int filecount, const char** filelist, uLong* crcvals, MPI_Comm comm)
{
  int rc = 0;

  int rank;
  MPI_Comm_rank(comm, &rank);

  double start = MPI_Wtime();

  redset d;
  int redset_rc = redset_recover(comm, path, &d);
  if (redset_rc != REDSET_SUCCESS) {
    ABORT("recover failed");
  }

  double end = MPI_Wtime();
  if (rank == 0) {
    printf("mode = %d: Recover time (lost 0 ranks): %f\n", mode, (end - start));
  }

  check_for_redundancy_files(mode, path, d);

  if (redset_delete(&d) != REDSET_SUCCESS) {
    ABORT("failed to delete redundancy descriptor");
  }

  check_crcs(filecount, filelist, crcvals);

  return rc;
}

int increment_index(int* index, int i, int k, int ranks)
{
  /* compute limit of current index */
  int limit = ranks - (k - i);
  if (index[i] < limit) {
    /* current index is below its limit, bump it up */
    index[i]++;
    return 0;
  } else if (i > 0) {
    /* current index has reached its limit,
     * try to bump the next index down if there is one */
    increment_index(index, i - 1, k, ranks);

    /* reset current index to one more than index below */
    index[i] = index[i - 1] + 1;
  } else {
    /* we've exhausted the lowest index value,
     * bump it up anyway so higher index values
     * also exceed their limits */
    index[i]++;
  }

  if (index[i] > limit) {
    /* we've exhausted the range for this value */
    return 1;
  }

  /* otherwise, we found a valid index value */
  return 0;
}

/* test that recover succeeds when no files are lost,
 * nothing to do, but should still return success */
int test_recover_loss_k_ranks(int mode, int apply_k, const char* path, int count, const char** filelist, uLong* crcvals, int k, MPI_Comm comm)
{
  int rc = 0;
  int tmp_rc;

  int rank, ranks;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &ranks);

  int* index = (int*) malloc(k * sizeof(int));

  /* initialize index values to 0, 1, ..., k-1 */
  int i;
  for (i = 0; i < k; i++) {
    index[i] = i;
  }

  int done = 0;
  while (! done) {
    /* delete files on rank 0 */
    if (rank == 0) {
      printf("mode = %d: deleting ranks = ", mode); 
      for (i = 0; i < k; i++) {
        printf("%d, ", index[i]); 
      }
      printf("\n");
      fflush(stdout);
    }
    for (i = 0; i < k; i++) {
      if (rank == index[i]) {
        delete_files(count, filelist);
      }
    }
  
    double start = MPI_Wtime();

    redset d;
    int redset_rc = redset_recover(comm, path, &d);
    int recovered = alltrue(redset_rc == REDSET_SUCCESS, comm);
  
    double end = MPI_Wtime();
    if (rank == 0) {
      printf("mode = %d: Recover time (protect %d ranks, lost %d ranks): %f\n", mode, apply_k, k, (end - start));
    }

    if ((mode == REDSET_COPY_SINGLE  && k > 0)       ||
        (mode == REDSET_COPY_PARTNER && k > apply_k) ||
        (mode == REDSET_COPY_XOR     && k > 1)       ||
        (mode == REDSET_COPY_RS      && k > apply_k))
     {
      if (recovered) {
        ABORT("should not have be able to recover files");
      }
    } else {
      if (recovered) {
        check_crcs(count, filelist, crcvals);
      } else {
        ABORT("recover failed");
      }
  
      check_for_redundancy_files(mode, path, d);
    }
  
    if (rc == 0 && recovered) {
      /* delete files on rank 0 */
      if (rank == 0) {
        printf("mode = %d: deleting ranks and redundancy = ", mode); 
        for (i = 0; i < k; i++) {
          printf("%d, ", index[i]); 
        }
        printf("\n");
        fflush(stdout);
      }
      for (i = 0; i < k; i++) {
        if (rank == index[i]) {
          delete_files(count, filelist);
          delete_redundancy_files(mode, path, d);
        }
      }
    }

    if (redset_delete(&d) != REDSET_SUCCESS) {
      ABORT("failed to delete redundancy descriptor");
    }

    /* current rank failed to recover, so no need to test other ranks */
    if (rc || !recovered) {
      fflush(stdout);
      return 1;
    }

    redset_rc = redset_recover(comm, path, &d);
    recovered = alltrue(redset_rc == REDSET_SUCCESS, comm);
  
    if ((mode == REDSET_COPY_SINGLE  && k > 0)       ||
        (mode == REDSET_COPY_PARTNER && k > apply_k) ||
        (mode == REDSET_COPY_XOR     && k > 1)       ||
        (mode == REDSET_COPY_RS      && k > apply_k))
    {
      if (recovered) {
        ABORT("should not have be able to recover files");
      }
    } else {
      if (recovered) {
        check_crcs(count, filelist, crcvals);
      } else {
        ABORT("recover failed");
      }
  
      check_for_redundancy_files(mode, path, d);
    }
  
    if (redset_delete(&d) != REDSET_SUCCESS) {
      ABORT("failed to delete redundancy descriptor");
    }

    /* current rank failed to recover, so no need to test other ranks */
    if (rc || !recovered) {
      fflush(stdout);
      return 1;
    }

    done = increment_index(index, k-1, k, ranks);
  }

  free(index);

  return rc;
}

int test_sequence(int copymode, const char* group, int filecount, const char** filelist, const char* prefix)
{
  int rc = 0;

  uLong* crcvals = (uLong*) malloc(filecount * sizeof(uLong));

  int ranks;
  MPI_Comm_size(MPI_COMM_WORLD, &ranks);

  int protect_k;
  for (protect_k = 0; protect_k < ranks; protect_k++) {
    /* TODO: check that redset returns an error on create for others */
    /* only single can handle k=0 protection */
    if (protect_k == 0 && copymode != REDSET_COPY_SINGLE) {
      continue;
    }

    /* no need to bother copy modes that do not support level of protection */
    if ((protect_k > 0 && copymode == REDSET_COPY_SINGLE) ||
        (protect_k > 1 && copymode == REDSET_COPY_XOR))
    {
      continue;
    }

    int lose_k;
    for (lose_k = 0; lose_k < ranks; lose_k++) {
      create_files(filecount, filelist);
      init_crcs(filecount, filelist, crcvals);

      redset d;
      switch (copymode) {
      case REDSET_COPY_SINGLE:
        redset_create_single(MPI_COMM_WORLD, group, &d);
        break;
      case REDSET_COPY_PARTNER:
        redset_create_partner(MPI_COMM_WORLD, group, 8, protect_k, &d);
        break;
      case REDSET_COPY_XOR:
        redset_create_xor(MPI_COMM_WORLD, group, 8, &d);
        break;
      case REDSET_COPY_RS:
        redset_create_rs(MPI_COMM_WORLD, group, 8, protect_k, &d);
        break;
      }

      test_apply(copymode, protect_k, filecount, filelist, crcvals, prefix, d, MPI_COMM_WORLD);

      if (lose_k == 0) {
        test_recover_no_loss(copymode, prefix, filecount, filelist, crcvals, MPI_COMM_WORLD);
      } else {
        test_recover_loss_k_ranks(copymode, protect_k, prefix, filecount, filelist, crcvals, lose_k, MPI_COMM_WORLD);
      }

      test_unapply(copymode, prefix, d, MPI_COMM_WORLD);

      int tmp_rc = redset_delete(&d);
      if (tmp_rc != REDSET_SUCCESS) {
        ABORT("failed to delete redundancy scheme object");
      }

      delete_files(filecount, filelist);
    }
  }

  free(crcvals);

  return rc;
}

int main (int argc, char* argv[])
{
  int tmp_rc;

  MPI_Init(&argc, &argv);

  int rank, ranks;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &ranks);

  char filename[256];
  sprintf(filename, "./testfile_%d.out", rank);

  int filecount = 1;
  const char* filelist[1] = { filename };

  char path[] = "scratch";
  mkdir(path, S_IRWXU);

  char prefix[256];
  sprintf(prefix, "%s/%d", path, rank);

  char hostname[1024];
  gethostname(hostname, sizeof(hostname));

  tmp_rc = redset_init();
  if (tmp_rc != REDSET_SUCCESS) {
    ABORT("redset_init failed");
  }

  test_sequence(REDSET_COPY_SINGLE,  hostname, filecount, filelist, prefix);
  test_sequence(REDSET_COPY_PARTNER, hostname, filecount, filelist, prefix);
  test_sequence(REDSET_COPY_XOR,     hostname, filecount, filelist, prefix);
  test_sequence(REDSET_COPY_RS,      hostname, filecount, filelist, prefix);

  tmp_rc = redset_finalize();
  if (tmp_rc != REDSET_SUCCESS) {
    ABORT("redset_init failed");
  }

  rmdir(path);

  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0) {
    printf("Done\n");
    fflush(stdout);
  }

  MPI_Finalize();

  return 0;
}
