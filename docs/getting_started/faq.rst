.. currentmodule:: hydromt_delft3dfm

.. _faq:

Frequently asked questions
==========================

This page contains some FAQ / tips and tricks to work with HydroMT-Delft3DFM.
For more general questions on how to work with data or the HydroMT config and command line,
you can visit the `HydroMT core FAQ page <https://deltares.github.io/hydromt/latest/getting_started/faq.html>`_

Building a dflowfm model
------------------------


Updating a dflowfm model
------------------------

 | **Q**: Can I select a specific dflowfm MDU config file when updating my model ?

It is possible. In that case, you need to start your HydroMT configuration with a **global** section
where you can specify which MDU file to use using the *config_fn* argument.

Others
------

