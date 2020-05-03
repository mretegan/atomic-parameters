The atomic structure codes originally developed by R. D. Cowan can be used to calculate the parameters of atomic configurations. The codes are however complicated to use. The current Python module addresses this by automatically creating the input files, running the programs, and extracting the parameters from the output files.

To use the script run the following commands:

.. code-block:: bash

    git clone https://github.com/mretegan/atomic-parameters.git
    cd atomic-parameters
    pip3 install --user -r requirements.txt
    python3 parameters.py --element "Fe" --configuration "1s1,3d5"

To run the tests go to the parent folder and run:

.. code-block:: bash

    python3 -m atomic-parameters.tests.test_parameters

The script uses the ``rcn``, ``rcn2``, ``ttrcg`` programs from the local ``cowan`` folder. These are available for macOS and Linux operating systems. See the README.rst for more details. If you want to use your own binaries, set the ``$TTMULT`` environment variable to the folder containing them.

