#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "redset.h"
#include "redset_util.h"
#include "redset_internal.h"

#include "kvtree.h"
#include "kvtree_util.h"

/* helper function to check for known options */
void check_known_options(const kvtree* configured_values,
                         const char* known_options[])
{
  /* report all unknown options (typos?) */
  const kvtree_elem* elem;
  for (elem = kvtree_elem_first(configured_values);
       elem != NULL;
       elem = kvtree_elem_next(elem))
  {
    const char* key = kvtree_elem_key(elem);

    /* must be only one level deep, ie plain kev = value */
    const kvtree* elem_hash = kvtree_elem_hash(elem);
    if (kvtree_size(elem_hash) != 1) {
      printf("Element %s has unexpected number of values: %d", key,
           kvtree_size(elem_hash));
      exit(EXIT_FAILURE);
    }

    const kvtree* kvtree_first_elem_hash =
      kvtree_elem_hash(kvtree_elem_first(elem_hash));
    if (kvtree_size(kvtree_first_elem_hash) != 0) {
      printf("Element %s is not a pure value", key);
      exit(EXIT_FAILURE);
    }

    /* check against known options */
    const char** opt;
    int found = 0;
    for (opt = known_options; *opt != NULL; opt++) {
      if (strcmp(*opt, key) == 0) {
        found = 1;
        break;
      }
    }
    if (! found) {
      printf("Unknown configuration parameter '%s' with value '%s'\n",
        kvtree_elem_key(elem),
        kvtree_elem_key(kvtree_elem_first(kvtree_elem_hash(elem)))
      );
      exit(EXIT_FAILURE);
    }
  }
}

/* helper function to check option values in kvtree against expected values */
void check_options(const int exp_debug, const int exp_set_size, const int exp_mpi_buf_size)
{
  kvtree* config = redset_config(NULL);
  if (config == NULL) {
    printf("redset_config failed\n");
    exit(EXIT_FAILURE);
  }

  int cfg_debug;
  if (kvtree_util_get_int(config, REDSET_KEY_CONFIG_DEBUG, &cfg_debug) !=
    KVTREE_SUCCESS)
  {
    printf("Could not get %s from redset_config\n",
           REDSET_KEY_CONFIG_DEBUG);
    exit(EXIT_FAILURE);
  }
  if (cfg_debug != exp_debug) {
    printf("redset_config returned unexpected value %d for %s. Expected %d.\n",
           cfg_debug, REDSET_KEY_CONFIG_DEBUG,
           exp_debug);
    exit(EXIT_FAILURE);
  }

  int cfg_set_size;
  if (kvtree_util_get_int(config, REDSET_KEY_CONFIG_SET_SIZE, &cfg_set_size) !=
    KVTREE_SUCCESS)
  {
    printf("Could not get %s from redset_config\n",
           REDSET_KEY_CONFIG_SET_SIZE);
    exit(EXIT_FAILURE);
  }
  if (cfg_set_size != exp_set_size) {
    printf("redset_config returned unexpected value %d for %s. Expected %d.\n",
           cfg_set_size, REDSET_KEY_CONFIG_SET_SIZE,
           exp_set_size);
    exit(EXIT_FAILURE);
  }

  int cfg_mpi_buf_size;
  if (kvtree_util_get_int(config, REDSET_KEY_CONFIG_MPI_BUF_SIZE,
                          &cfg_mpi_buf_size) != KVTREE_SUCCESS)
  {
    printf("Could not get %s from redset_config\n",
           REDSET_KEY_CONFIG_MPI_BUF_SIZE);
    exit(EXIT_FAILURE);
  }
  if (cfg_mpi_buf_size != exp_mpi_buf_size) {
    printf("redset_config returned unexpected value %d for %s. Expected %d.\n",
           cfg_mpi_buf_size, REDSET_KEY_CONFIG_MPI_BUF_SIZE,
           exp_mpi_buf_size);
    exit(EXIT_FAILURE);
  }

  static const char* known_options[] = {
    REDSET_KEY_CONFIG_DEBUG,
    REDSET_KEY_CONFIG_SET_SIZE,
    REDSET_KEY_CONFIG_MPI_BUF_SIZE,
    NULL
  };
  check_known_options(config, known_options);

  kvtree_delete(&config);
}

int
main(int argc, char *argv[]) {
    int rc;
    kvtree* redset_config_values = kvtree_new();

    MPI_Init(&argc, &argv);

    rc = redset_init();
    if (rc != REDSET_SUCCESS) {
        printf("redset_init() failed (error %d)\n", rc);
        return rc;
    }

    int old_redset_debug = redset_debug;
    int old_redset_set_size = redset_set_size;
    int old_redset_mpi_buf_size = redset_mpi_buf_size;

    int new_redset_debug = !old_redset_debug;
    int new_redset_set_size = old_redset_set_size + 1;
    int new_redset_mpi_buf_size = old_redset_mpi_buf_size + 1;

    check_options(old_redset_debug, old_redset_set_size,
                  old_redset_mpi_buf_size);

    /* check redset configuration settings */
    rc = kvtree_util_set_int(redset_config_values, REDSET_KEY_CONFIG_DEBUG,
                             new_redset_debug);
    if (rc != KVTREE_SUCCESS) {
        printf("kvtree_util_set_int failed (error %d)\n", rc);
        return rc;
    }
    rc = kvtree_util_set_int(redset_config_values, REDSET_KEY_CONFIG_SET_SIZE,
                             new_redset_set_size);
    if (rc != KVTREE_SUCCESS) {
        printf("kvtree_util_set_int failed (error %d)\n", rc);
        return rc;
    }

    printf("Configuring redset (first set of options)...\n");
    if (redset_config(redset_config_values) == NULL) {
        printf("redset_config() failed\n");
        return EXIT_FAILURE;
    }

    if (redset_debug != new_redset_debug) {
        printf("redset_config() failed to set %s: %d != %d\n",
               REDSET_KEY_CONFIG_DEBUG, redset_debug, new_redset_debug);
        return EXIT_FAILURE;
    }

    if (redset_set_size != new_redset_set_size) {
        printf("REDSET_Config() failed to set %s: %d != %d\n",
               REDSET_KEY_CONFIG_SET_SIZE, redset_set_size, old_redset_set_size);
        return EXIT_FAILURE;
    }

    check_options(new_redset_debug, new_redset_set_size,
                  old_redset_mpi_buf_size);

    /* configure remaining options */
    kvtree_delete(&redset_config_values);
    redset_config_values = kvtree_new();

    rc = kvtree_util_set_int(redset_config_values, REDSET_KEY_CONFIG_MPI_BUF_SIZE,
                             new_redset_mpi_buf_size);
    if (rc != KVTREE_SUCCESS) {
        printf("kvtree_util_set_int failed (error %d)\n", rc);
        return rc;
    }

    printf("Configuring redset (second set of options)...\n");
    if (redset_config(redset_config_values) == NULL) {
        printf("redset_config() failed\n");
        return EXIT_FAILURE;
    }

    if (redset_debug != new_redset_debug) {
        printf("redset_config() failed to set %s: %d != %d\n",
               REDSET_KEY_CONFIG_DEBUG, redset_debug, new_redset_debug);
        return EXIT_FAILURE;
    }

    if (redset_set_size != new_redset_set_size) {
        printf("REDSET_Config() failed to set %s: %d != %d\n",
               REDSET_KEY_CONFIG_SET_SIZE, redset_set_size, old_redset_set_size);
        return EXIT_FAILURE;
    }

    if (redset_mpi_buf_size != new_redset_mpi_buf_size) {
        printf("REDSET_Config() failed to set %s: %d != %d\n",
               REDSET_KEY_CONFIG_MPI_BUF_SIZE, redset_mpi_buf_size,
               old_redset_mpi_buf_size);
        return EXIT_FAILURE;
    }

    check_options(new_redset_debug, new_redset_set_size,
                  new_redset_mpi_buf_size);

    rc = redset_finalize();
    if (rc != REDSET_SUCCESS) {
        printf("redset_finalize() failed (error %d)\n", rc);
        return rc;
    }

    MPI_Finalize();

    return REDSET_SUCCESS;
}
