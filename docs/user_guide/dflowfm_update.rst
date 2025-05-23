.. _delft3dfm_update:

Updating a model
----------------
To add or change one or more components of an existing Delft3D FM model the ``update`` method can be used.

**Steps in brief:**

1) You have an **existing model** schematization. This model does not have to be complete.
2) Prepare or use a pre-defined **data catalog** with all the required data sources, see `working with data <https://deltares.github.io/hydromt/latest/user_guide/data_prepare_cat.html>`_.
3) Prepare a **model configuration** with the methods that you want to use to add or change components of your model: see `model configuration <https://deltares.github.io/hydromt/latest/user_guide/model_config.html>`_.
4) **Update** your model using the CLI or Python interface.

.. code-block:: console

    activate hydromt-delft3dfm
    hydromt update delft3dfm path/to/model_to_update -o path/to/updated_model -i delft3dfm_update.yml -d data_sources.yml -vvv

.. NOTE::

    By default, the updated model will overwrite your existing one. To save the updated model in a different
    folder, use the -o path/to/updated_model option of the CLI.

.. TIP::

    By default all model data is written at the end of the update method. If your update however
    only affects a certain model data (e.g. staticmaps or forcing) you can add a write_* method
    (e.g. `write_maps`, `write_forcing`) to the .yml file and only these data will be written.

    Note that the model config is often changed as part of the a model method and `write_config`
    should thus often be added to the .yml file to keep the model data and config consistent.

.. toctree::
    :hidden:

    Example: Update Delft3D FM model (refine 2dgrid) <../_examples/update_refine_2dgrid.ipynb>
