#include <stdlib.h>
#include <stdio.h>
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

/* compute MPI_LAND of value across all procs,
 * return 1 if value is 1 on all procs, 0 otherwise */
int alltrue(int value, MPI_Comm comm)
{
  int all_value;
  MPI_Allreduce(&value, &all_value, 1, MPI_INT, MPI_LAND, comm);
  return all_value;
}

void create_files(int count, const char** filelist)
{
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
      printf("Error opening file %s: %d %s\n", name, errno, strerror(errno));
    }
  }

  free(buf);
  return;
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
      printf("ERROR: Failed CRC32 check for %s\n", file);
      return 1;
    }
  }
  return 0;
}

/* delete the list of files */
void delete_files(int count, const char** filelist)
{
  int i;
  for (i = 0; i < count; i++) {
    const char* name = filelist[i];
    unlink(name);
  }
  return;
}

/* check that we have each redundancy file */
int check_for_redundancy_files(int mode, const char* path, redset d)
{
  int rc = 0;

  /* get list of redundancy files */
  redset_filelist list = redset_filelist_get(path, d);
  if (list == NULL) {
    printf("ERROR: Failed to get list of redundancy files\n");
    return 1;
  }

  /* check that each file in the list exists */
  int i;
  int count = redset_filelist_count(list);
  for (i = 0; i < count; i++) {
    const char* name = redset_filelist_file(list, i);
    if (access(name, F_OK) != 0) {
      printf("ERROR: Missing redundancy file %s\n", name);
      rc = 1;
    }
  }

  /* check that we got the expected number of files in the list */
  if (count != 2) {
    /* all methods generate two files */
    printf("ERROR: Unexpected number of redundancy files\n");
    rc = 1;
  }

  /* free the list */
  redset_filelist_release(&list);
  if (list != NULL) {
    printf("ERROR: Failed to free list of redundancy files\n");
    rc = 1;
  }

  return rc;
}

int delete_redundancy_files(int mode, const char* path, redset d)
{
  int rc = 0;

  /* get list of redundancy files */
  redset_filelist list = redset_filelist_get(path, d);
  if (list == NULL) {
    printf("ERROR: Failed to get list of redundancy files\n");
    return 1;
  }

  /* check that each file in the list exists */
  int i;
  int count = redset_filelist_count(list);
  for (i = 0; i < count; i++) {
    const char* name = redset_filelist_file(list, i);
    if (unlink(name) != 0) {
      printf("ERROR: Deleting redundancy file %s\n", name);
      rc = 1;
    }
  }

  /* check that we got the expected number of files in the list */
  if (count != 2) {
    /* all methods generate two files */
    printf("ERROR: Unexpected number of redundancy files\n");
    rc = 1;
  }

  /* free the list */
  redset_filelist_release(&list);
  if (list != NULL) {
    printf("ERROR: Failed to free list of redundancy files\n");
    rc = 1;
  }

  return rc;
}

/* apply redundancy descriptor and check that redundancy files are created */
int test_apply(int mode, int filecount, const char** filelist, const char* path, redset d, MPI_Comm comm)
{
  int rc = 0;

  int rank;
  MPI_Comm_rank(comm, &rank);

  double start = MPI_Wtime();

  int redset_rc = redset_apply(filecount, filelist, path, d);
  if (! alltrue(redset_rc == REDSET_SUCCESS, comm)) {
    printf("ERROR: apply failed\n");
    rc = 1;
  }

  double end = MPI_Wtime();
  if (rank == 0) {
    printf("mode = %d: Apply time: %f\n", mode, (end - start));
  }

  rc = check_for_redundancy_files(mode, path, d);

  return rc;
}

/* unapply redundancy descriptor and check that redundancy files are gone */
int test_unapply(int mode, const char* path, redset d, MPI_Comm comm)
{
  int rc = 0;

  int redset_rc = redset_unapply(path, d);
  if (! alltrue(redset_rc == REDSET_SUCCESS, comm)) {
    printf("ERROR: unapply failed\n");
    rc = 1;
  }

  /* get list of redundancy files */
  redset_filelist list = redset_filelist_get(path, d);
  if (list == NULL) {
    printf("ERROR: Failed to get list of redundancy files\n");
    return 1;
  }

  /* check that each file in the list is gone */
  int i;
  int count = redset_filelist_count(list);
  for (i = 0; i < count; i++) {
    const char* name = redset_filelist_file(list, i);
    if (access(name, F_OK) == 0) {
      printf("ERROR: Found redundancy file %s\n", name);
      rc = 1;
    }
  }

  /* check that we got the expected number of files in the list */
  if (count != 2) {
    /* all methods generate two files */
    printf("ERROR: Unexpected number of redundancy files\n");
    rc = 1;
  }

  /* free the list */
  redset_filelist_release(&list);
  if (list != NULL) {
    printf("ERROR: Failed to free list of redundancy files\n");
    rc = 1;
  }

  return rc;
}

/* test that recover succeeds when no files are lost,
 * nothing to do, but should still return success */
int test_recover_no_loss(int mode, const char* path, int filecount, const char** filelist, MPI_Comm comm)
{
  int rc = 0;

  int rank;
  MPI_Comm_rank(comm, &rank);

  double start = MPI_Wtime();

  redset d;
  int redset_rc = redset_recover(comm, path, &d);
  if (! alltrue(redset_rc == REDSET_SUCCESS, comm)) {
    printf("ERROR: recover failed\n");
    rc = 1;
  }

  double end = MPI_Wtime();
  if (rank == 0) {
    printf("mode = %d: Recover time (no loss): %f\n", mode, (end - start));
  }

  rc = check_for_redundancy_files(mode, path, d);

  if (redset_delete(&d) != REDSET_SUCCESS) {
    printf("ERROR: failed to delete redundancy descriptor\n");
    rc = 1;
  }

  return rc;
}

/* test that recover succeeds when no files are lost,
 * nothing to do, but should still return success */
int test_recover_loss_one_rank(int mode, const char* path, int count, const char** filelist, uLong* crcvals, MPI_Comm comm)
{
  int rc = 0;

  int rank, ranks;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &ranks);

  int r;
  for (r = 0; r < ranks; r++) {
    /* delete files on rank 0 */
    if (rank == 0) {
      printf("mode = %d: deleting rank = %d\n", mode, r);  fflush(stdout);
    }
    if (rank == r) {
      delete_files(count, filelist);
    }

    double start = MPI_Wtime();

    redset d;
    int redset_rc = redset_recover(comm, path, &d);
    int recovered = alltrue(redset_rc == REDSET_SUCCESS, comm);

    double end = MPI_Wtime();
    if (rank == 0) {
      printf("mode = %d: Recover time (one rank): %f\n", mode, (end - start));
    }

    if (mode == REDSET_COPY_SINGLE) {
      if (recovered) {
        printf("ERROR: should not have be able to recover files\n");
        rc = 1;
      }
    } else {
      if (recovered) {
        rc = check_crcs(count, filelist, crcvals);
      } else {
        printf("ERROR: recover failed\n");
        rc = 1;
      }

      int tmp_rc = check_for_redundancy_files(mode, path, d);
      if (tmp_rc) {
        rc = tmp_rc;
      }
    }

    if (rc == 0 && recovered) {
      if (rank == 0) {
        printf("mode = %d: deleting rank = %d and its redundancy files\n", mode, r);  fflush(stdout);
      }
      if (rank == r) {
        delete_files(count, filelist);
        delete_redundancy_files(mode, path, d);
      }
    }

    if (redset_delete(&d) != REDSET_SUCCESS) {
      printf("ERROR: failed to delete redundancy descriptor\n");
      rc = 1;
    }

    /* current rank failed to recover, so no need to test other ranks */
    if (rc || !recovered) {
      fflush(stdout);
      return 1;
    }

    redset_rc = redset_recover(comm, path, &d);
    recovered = alltrue(redset_rc == REDSET_SUCCESS, comm);

    if (mode == REDSET_COPY_SINGLE) {
      if (recovered) {
        printf("ERROR: should not have be able to recover files\n");
        rc = 1;
      }
    } else {
      if (recovered) {
        rc = check_crcs(count, filelist, crcvals);
      } else {
        printf("ERROR: recover failed\n");
        rc = 1;
      }

      int tmp_rc = check_for_redundancy_files(mode, path, d);
      if (tmp_rc) {
        rc = tmp_rc;
      }
    }

    if (redset_delete(&d) != REDSET_SUCCESS) {
      printf("ERROR: failed to delete redundancy descriptor\n");
      rc = 1;
    }

    /* current rank failed to recover, so no need to test other ranks */
    if (rc || !recovered) {
      fflush(stdout);
      return 1;
    }
  }

  return rc;
}

/* test that recover succeeds when no files are lost,
 * nothing to do, but should still return success */
int test_recover_loss_two_ranks(int mode, const char* path, int count, const char** filelist, uLong* crcvals, MPI_Comm comm)
{
  int rc = 0;

  int rank, ranks;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &ranks);

  int i, j;
  for (i = 0; i < ranks; i++) {
  for (j = i + 1; j < ranks; j++) {
    /* delete files on rank 0 */
    if (rank == 0) {
      printf("mode = %d: deleting ranks = %d,%d\n", mode, i, j);  fflush(stdout);
    }
    if (rank == i || rank == j) {
      delete_files(count, filelist);
    }
  
    double start = MPI_Wtime();

    redset d;
    int redset_rc = redset_recover(comm, path, &d);
    int recovered = alltrue(redset_rc == REDSET_SUCCESS, comm);
  
    double end = MPI_Wtime();
    if (rank == 0) {
      printf("mode = %d: Recover time (two ranks): %f\n", mode, (end - start));
    }

    if (mode == REDSET_COPY_SINGLE || mode == REDSET_COPY_XOR) {
      if (recovered) {
        printf("ERROR: should not have be able to recover files\n");
        rc = 1;
      }
    } else {
      if (recovered) {
        rc = check_crcs(count, filelist, crcvals);
      } else {
        printf("ERROR: recover failed\n");
        rc = 1;
      }
  
      int tmp_rc = check_for_redundancy_files(mode, path, d);
      if (tmp_rc) {
        rc = tmp_rc;
      }
    }
  
    if (rc == 0 && recovered) {
      if (rank == 0) {
        printf("mode = %d: deleting ranks = %d,%d and their redundancy files\n", mode, i, j);  fflush(stdout);
      }
      if (rank == i || rank == j) {
        delete_files(count, filelist);
        delete_redundancy_files(mode, path, d);
      }
    }

    if (redset_delete(&d) != REDSET_SUCCESS) {
      printf("ERROR: failed to delete redundancy descriptor\n");
      rc = 1;
    }

    /* current rank failed to recover, so no need to test other ranks */
    if (rc || !recovered) {
      fflush(stdout);
      return 1;
    }

    redset_rc = redset_recover(comm, path, &d);
    recovered = alltrue(redset_rc == REDSET_SUCCESS, comm);
  
    if (mode == REDSET_COPY_SINGLE || mode == REDSET_COPY_XOR) {
      if (recovered) {
        printf("ERROR: should not have be able to recover files\n");
        rc = 1;
      }
    } else {
      if (recovered) {
        rc = check_crcs(count, filelist, crcvals);
      } else {
        printf("ERROR: recover failed\n");
        rc = 1;
      }
  
      int tmp_rc = check_for_redundancy_files(mode, path, d);
      if (tmp_rc) {
        rc = tmp_rc;
      }
    }
  
    if (redset_delete(&d) != REDSET_SUCCESS) {
      printf("ERROR: failed to delete redundancy descriptor\n");
      rc = 1;
    }

    /* current rank failed to recover, so no need to test other ranks */
    if (rc || !recovered) {
      fflush(stdout);
      return 1;
    }
  }
  }

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
int test_recover_loss_k_ranks(int mode, const char* path, int count, const char** filelist, uLong* crcvals, int k, MPI_Comm comm)
{
  int rc = 0;

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
      printf("mode = %d: Recover time (two ranks): %f\n", mode, (end - start));
    }

    if ((mode == REDSET_COPY_SINGLE && k > 0) ||
        (mode == REDSET_COPY_XOR    && k > 1))
     {
      if (recovered) {
        printf("ERROR: should not have be able to recover files\n");
        rc = 1;
      }
    } else {
      if (recovered) {
        rc = check_crcs(count, filelist, crcvals);
      } else {
        printf("ERROR: recover failed\n");
        rc = 1;
      }
  
      int tmp_rc = check_for_redundancy_files(mode, path, d);
      if (tmp_rc) {
        rc = tmp_rc;
      }
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
      printf("ERROR: failed to delete redundancy descriptor\n");
      rc = 1;
    }

    /* current rank failed to recover, so no need to test other ranks */
    if (rc || !recovered) {
      fflush(stdout);
      return 1;
    }

    redset_rc = redset_recover(comm, path, &d);
    recovered = alltrue(redset_rc == REDSET_SUCCESS, comm);
  
    if ((mode == REDSET_COPY_SINGLE && k > 0) ||
        (mode == REDSET_COPY_XOR    && k > 1))
    {
      if (recovered) {
        printf("ERROR: should not have be able to recover files\n");
        rc = 1;
      }
    } else {
      if (recovered) {
        rc = check_crcs(count, filelist, crcvals);
      } else {
        printf("ERROR: recover failed\n");
        rc = 1;
      }
  
      int tmp_rc = check_for_redundancy_files(mode, path, d);
      if (tmp_rc) {
        rc = tmp_rc;
      }
    }
  
    if (redset_delete(&d) != REDSET_SUCCESS) {
      printf("ERROR: failed to delete redundancy descriptor\n");
      rc = 1;
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

void test_sequence(int copymode, const char* group, int filecount, const char** filelist, const char* prefix)
{
  create_files(filecount, filelist);

  uLong* crcvals = (uLong*) malloc(filecount * sizeof(uLong));
  init_crcs(filecount, filelist, crcvals);

  redset d;
  //redset_create(copymode, MPI_COMM_WORLD, group, &d);
  switch (copymode) {
  case REDSET_COPY_SINGLE:
    redset_create_single(MPI_COMM_WORLD, group, &d);
    break;
  case REDSET_COPY_PARTNER:
    redset_create_partner(MPI_COMM_WORLD, group, 8, 1, &d);
    break;
  case REDSET_COPY_XOR:
    redset_create_xor(MPI_COMM_WORLD, group, 8, &d);
    break;
  case REDSET_COPY_RS:
    redset_create_rs(MPI_COMM_WORLD, group, 8, 2, &d);
    break;
  }

  test_apply(copymode, filecount, filelist, prefix, d, MPI_COMM_WORLD);
  check_crcs(filecount, filelist, crcvals);

  test_recover_no_loss(copymode, prefix, filecount, filelist, MPI_COMM_WORLD);
  check_crcs(filecount, filelist, crcvals);

  test_recover_loss_one_rank(copymode, prefix, filecount, filelist, crcvals, MPI_COMM_WORLD);

  test_unapply(copymode, prefix, d, MPI_COMM_WORLD);

  redset_delete(&d);

  delete_files(filecount, filelist);



  int ranks;
  MPI_Comm_size(MPI_COMM_WORLD, &ranks);

  int k;
  for (k = 2; k < ranks; k++) {
    create_files(filecount, filelist);
    init_crcs(filecount, filelist, crcvals);

    //redset_create(copymode, MPI_COMM_WORLD, group, &d);
    switch (copymode) {
    case REDSET_COPY_SINGLE:
      redset_create_single(MPI_COMM_WORLD, group, &d);
      break;
    case REDSET_COPY_PARTNER:
      redset_create_partner(MPI_COMM_WORLD, group, 8, k, &d);
      break;
    case REDSET_COPY_XOR:
      redset_create_xor(MPI_COMM_WORLD, group, 8, &d);
      break;
    case REDSET_COPY_RS:
      redset_create_rs(MPI_COMM_WORLD, group, 8, k, &d);
      break;
    }

    test_apply(copymode, filecount, filelist, prefix, d, MPI_COMM_WORLD);
    check_crcs(filecount, filelist, crcvals);

    //test_recover_loss_two_ranks(copymode, prefix, filecount, filelist, crcvals, MPI_COMM_WORLD);
    test_recover_loss_k_ranks(copymode, prefix, filecount, filelist, crcvals, k, MPI_COMM_WORLD);

    test_unapply(copymode, prefix, d, MPI_COMM_WORLD);

    redset_delete(&d);

    delete_files(filecount, filelist);
  }

  free(crcvals);

  return;
}

int main (int argc, char* argv[])
{
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

  redset_init();

  char hostname[1024];
  gethostname(hostname, sizeof(hostname));

  test_sequence(REDSET_COPY_SINGLE,  hostname, filecount, filelist, prefix);
  test_sequence(REDSET_COPY_PARTNER, hostname, filecount, filelist, prefix);
  test_sequence(REDSET_COPY_XOR,     hostname, filecount, filelist, prefix);
  test_sequence(REDSET_COPY_RS,      hostname, filecount, filelist, prefix);

  redset_finalize();

  rmdir(path);

  MPI_Finalize();

  return 0;
}
