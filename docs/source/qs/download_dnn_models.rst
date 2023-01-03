Download DNN Models
======================================
The neural network models used in the tutorial examples are indepentently trained
by our collaborators team – `DeepCombustion <https://github.com/deepcombustion/deepcombustion>`_. 
To run DeepFlame with DNN, first download the DeepCombustion repository into ``deepflame-dev/``: 

.. code-block:: bash

    cd $DF_ROOT
    git clone https://github.com/deepcombustion/deepcombustion.git 
    
Then copy the required DNN model into ``mechanisms/``, for example:

.. code-block:: bash

    cp -r deepcombustion/DeePCK/Model/HE04_Hydrogen_ESH2_GMS_sub_20221101/ mechanisms/
    
.. Note:: Here ``HE04_Hydrogen_ESH2_GMS_sub_20221101`` is the default DNN model for all the tutorial cases in ``$DF_ROOT/examples/``.
