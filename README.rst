The atomic structure codes originally developed by R. D. Cowan can be used to calculate the parameters of atomic configurations. The codes are however complicated to use. The current Python module addresses this by automatically creating the input files, running the programs, and extracting the parameters from the output files.

To use the script run the following commands:

.. code-block:: bash

    git clone https://github.com/mretegan/atomic-parameters.git
    cd atomic-parameters
    python parameters.py --element "Fe" --configuration "1s1,3d5"

Note that the programs ``rcn``, ``rcn2``, ``ttrcg`` must be installed and the ``$TTMULT_HOME`` environment variable must be set to the folder containing the programs.

