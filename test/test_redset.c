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

  int i;
  for (i = 0; i < count; i++) {
    char buf[256];
    sprintf(buf, "data from rank %d in file %d\n", rank, i);

    const char* name = filelist[i];
    int fd = open(name, O_WRONLY | O_TRUNC | O_CREAT, S_IRUSR | S_IWUSR);
    if (fd != -1) {
      write(fd, buf, strlen(buf));
      close(fd);
    } else {
      printf("Error opening file %s: %d %s\n", name, errno, strerror(errno));
    }
  }
  return;
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

/* apply redundancy descriptor and check that redundancy files are created */
int test_apply(int mode, int filecount, const char** filelist, const char* path, redset d, MPI_Comm comm)
{
  int rc = 0;

  int redset_rc = redset_apply(filecount, filelist, path, d);
  if (! alltrue(redset_rc == REDSET_SUCCESS, comm)) {
    printf("ERROR: apply failed\n");
    rc = 1;
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

  redset d;
  int redset_rc = redset_recover(comm, path, &d);
  if (! alltrue(redset_rc == REDSET_SUCCESS, comm)) {
    printf("ERROR: recover failed\n");
    rc = 1;
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
int test_recover_loss_one_rank(int mode, const char* path, int count, const char** filelist, MPI_Comm comm)
{
  int rc = 0;

  int rank, ranks;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &ranks);

  /* delete files on rank 0 */
  if (rank == 0) {
    delete_files(count, filelist);
  }

  redset d;
  int redset_rc = redset_recover(comm, path, &d);
  int recovered = alltrue(redset_rc == REDSET_SUCCESS, comm);

  if (mode == REDSET_COPY_SINGLE) {
    if (recovered) {
      printf("ERROR: should not have be able to recover files\n");
      rc = 1;
    }
  } else {
    if (! recovered) {
      printf("ERROR: recover failed\n");
      rc = 1;
    }

    rc = check_for_redundancy_files(mode, path, d);
  }

  if (redset_delete(&d) != REDSET_SUCCESS) {
    printf("ERROR: failed to delete redundancy descriptor\n");
    rc = 1;
  }

  return rc;
}

void test_sequence(int copymode, const char* group, int filecount, const char** filelist, const char* prefix)
{
  create_files(filecount, filelist);

  redset d;
  redset_create(copymode, MPI_COMM_WORLD, group, &d);

  test_apply(copymode, filecount, filelist, prefix, d, MPI_COMM_WORLD);

  test_recover_no_loss(copymode, prefix, filecount, filelist, MPI_COMM_WORLD);

  test_recover_loss_one_rank(copymode, prefix, filecount, filelist, MPI_COMM_WORLD);

  test_unapply(copymode, prefix, d, MPI_COMM_WORLD);

  redset_delete(&d);

  delete_files(filecount, filelist);

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

  redset_finalize();

  rmdir(path);

  MPI_Finalize();

  return 0;
}
