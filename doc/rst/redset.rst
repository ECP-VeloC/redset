redset
======
A redundancy descriptor is a data structure that describes how a dataset is encoded.
It tracks information such as the list of source files and the redundancy scheme that is applied.
The data structure also records information on the group of processes that make up a redundancy set,
such as the number of processes in the set, as well as,
a unique integer that identifies the set, called the "group id".

This section describes some of the most common redundancy descriptor functions.

For more documentation see src/redset.h and the examples in test/test_redset.c.

For developers, some high-level implementation details are in doc/rst/implementation.rst.

Initializing and freeing redundancy descriptors
+++++++++++++++++++++++++++++++++++++++++++++++
Create a new redundancy set::

  redset d;
  redset_create(type, comm, group, &d);

The ``type`` specifies the redundancy scheme that should be applied.
It should be one of: ``REDSET_COPY_SINGLE``, ``REDSET_COPY_PARTNER``, ``REDSET_COPY_XOR``.

The ``comm`` value specifies the MPI communicator containing the processes participating in the set.

The ``group`` value is a string specifying the failure group.
All processes proving the same value for ``group`` are expected to fail at the same time.
For example, one might use the hostname to indicate that all processes on the same host are likely to fail simultaneously.

Free memory associated with a redundancy descriptor::

  redset_delete(&d);

Apply, recover, and unapply a redundancy scheme
+++++++++++++++++++++++++++++++++++++++++++++++
To apply a redundancy scheme to a set of files::

  redset_appy(numfiles, files, name, d);

Here, ``files`` should be an array of file names and ``numfiles`` provides the size of the array.
The ``name`` parameter specifies a prefix string to prepend to the name of all redundancy files the library creates.
The ``name`` string can include directory components to write redundancy files to a different directory,
otherwise they will be written to the current working directory.

To recover files accoding to a redundancy scheme::

  redset d;
  int ret = redset_recover(comm, name, &d);
  if (ret == REDSET_SUCCESS) {
    // data recovered
  } else {
    // data lost
  }

Here ``comm`` should correspond to the same parent communicator used to create the original redundancy descriptor,
and ``name`` should specify the prefix string used when the redundancy descriptor was applied.
If the rebuild succeeds, this will return REDSET_SUCCESS.
It returns the same value to all processes in ``comm``.

To remove a redundancy scheme and delete assocated redundancy files::

  redset_unapply(name, d);

Again, ``name`` should specify the string prefix used when the redundancy scheme was applied.

Listing redundancy files
++++++++++++++++++++++++
Sometimes it may be necessary to obtain a list of files added when the redundancy scheme was applied::

  int i;
  redset_filelist list = redset_filelist_get(name, d);
  int count = redset_filelist_count(list);
  for (i = 0; i < count; i++) {
    const char* file = redset_filelist_file(list, i);
  }
