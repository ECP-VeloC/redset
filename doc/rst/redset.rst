redset
======
A redundancy descriptor is a data structure that describes how a dataset should be encoded.
It tracks information including the redundancy scheme to be applied
and the group of processes that make up a redundancy set.

This section describes some of the most common redundancy descriptor functions.

For more documentation see src/redset.h and the examples in test/test_redset.c.

For developers, some high-level implementation details are in doc/rst/implementation.rst.

Create and delete redundancy descriptors
++++++++++++++++++++++++++++++++++++++++
Create a new redundancy descriptor::

  redset d;
  redset_create(type, comm, group, &d);

The ``type`` specifies the redundancy scheme that should be applied.
It should be one of: ``REDSET_COPY_SINGLE``, ``REDSET_COPY_PARTNER``, ``REDSET_COPY_XOR``, or ``REDSET_COPY_RS``.

The ``comm`` value specifies the MPI communicator containing the processes participating in the set.
All redset functions are implied collectives over the set of processes in ``comm``.

The ``group`` value is a string specifying the failure group.
All processes proving the same value for ``group`` are expected to fail at the same time.
For example, one might use the hostname to indicate that all processes on the same host are likely to fail simultaneously.

Alternatively, each redundancy scheme has a custom create method
to allow the caller to specify certain parameters specific to each scheme::

  redset_create_single(comm, group, &d);
  redset_create_partner(comm, group, set_size, replicas, &d);
  redset_create_xor(comm, group, set_size, &d);
  redset_create_rs(comm, group, set_size, k, &d);

The ``comm`` and ``group`` parameters are the same as described above.

The ``set_size`` parameter specifies the minimum number of processes to include in each redundancy set.
The redset library uses this parameter to divide large process groups into smaller sets
such that no subset has less than ``set_size`` processes.
This division is done after subsetting processes across failure groups.
If the original set of processes is less than ``set_size``, redset constructs the largest set possible.
To illustrate, if the caller specifies ``set_size=8``,
redset creates redundancy groups of the following sizes for each count of available processes::

  procs   resulting set sizes
  4       4    (create a set as large as possible when total procs are less than 8)
  8       8
  9       9    (cannot divide 9 so that each subset has at least 8 procs)
  15      15
  16      8, 8 (can divide 16 so that each subset has at least 8 procs)
  17      9, 8
  18      9, 9

The ``replicas`` parameter specifies the number of partner copies to make.
This can range from 1 to one less than the set size.

The ``k`` parameter specifies the number of encoding blocks for Reed-Solomon.
This can range from 1 to one less than the set size.
Each group can recover from up to ``k`` simultaneous failures within the redundancy set.

To free resources associated with a redundancy descriptor::

  redset_delete(&d);

Apply and unapply a redundancy scheme
+++++++++++++++++++++++++++++++++++++
To apply a redundancy scheme to a set of files::

  redset_apply(numfiles, files, name, d);

Here ``files`` is an array of file names and ``numfiles`` provides the size of the array.
The ``name`` parameter specifies a prefix string to prepend to the name of all redundancy files that redset creates.
The ``name`` string can include directory components to write redundancy files into a particular directory,
otherwise redundancy files are written to the current working directory.
When specifying a different directory, the target directory must already exist.

To remove a redundancy scheme and delete the associated redundancy files::

  redset_unapply(name, d);

The ``name`` parameter must specify the same string prefix used when the redundancy scheme was applied.

Recover files with a redundancy scheme
++++++++++++++++++++++++++++++++++++++
When recovering files, one provides communicator and name as input parameters,
and one receives a redundancy descriptor as output::

  redset d;
  int ret = redset_recover(comm, name, &d);
  if (ret == REDSET_SUCCESS) {
    // data recovered
  } else {
    // data lost
  }
  redset_delete(&d);

The ``comm`` parameter must be equivalent to the parent communicator that was used to create the original redundancy descriptor.
The ``name`` parameter must specify the same prefix string that was used when the redundancy descriptor was applied.
The call to ``redset_recover`` returns the same return code to all processes in ``comm``.
If the recover operation succeeds, ``redset_recover`` returns REDSET_SUCCESS.
The recover operation is successful if either redset verifies that all source and associated redundancy files still exist
or if redset can successfully rebuild any files that are missing.

Whether successful or not, ``redset_recover`` returns a redundancy descriptor in the output parameter ``d``.
One may call ``redset_unapply`` using the returned descriptor to attempt to remove redset redundancy files,
though if the recover operation fails, redset may not have sufficient information to remove all redundancy files.
The caller is responsible for freeing the descriptor returned in ``d`` with a call to ``redset_free``.

Listing redundancy files
++++++++++++++++++++++++
Sometimes it is necessary for the caller to obtain a list of files added when the redundancy scheme was applied::

  int i;
  redset_filelist list = redset_filelist_get(name, d);
  int count = redset_filelist_count(list);
  for (i = 0; i < count; i++) {
    const char* file = redset_filelist_file(list, i);
  }
  redset_filelist_release(&list);
