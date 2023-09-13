"""
A pert file template.
"""


def pert(t, system):
    """
    Perturbation function called at each step.

    The function needs to be named ``pert`` and takes two positional arguments:
    ``t`` for the simulation time, and ``system`` for the system object.
    Arbitrary logic and calculations can be applied in this function to
    ``system``.

    If the perturbation event involves switching, such as disconnecting a line,
    one will need to set the ``system.TDS.custom_event`` flag to ``True`` to
    trigger a system connectivity checking, and Jacobian rebuilding and
    refactorization. To implement, add the following line to the scope where the
    event is triggered:

    .. code-block :: python

        system.TDS.custom_event = True

    In other scopes of the code where events are not triggered, do not add the
    above line as it may cause significant slow-down.

    The perturbation file can be supplied to the CLI using the ``--pert``
    argument or supplied to :py:func:`andes.main.run` using the ``pert``
    keyword.

    Parameters
    ----------
    t : float
        Simulation time.
    system : andes.system.System
        System object supplied by the simulator.

    """

    pert_my(t, system)


def pert_my(t, system):
    """
    Perturbation function called at each step.

    Parameters
    ----------
    t : float
        Simulation time.
    system : andes.system.System
        System object supplied by the simulator.

    """

    # read the file only upon the first run

    if system.dae.t == 0:  
        # upon initialization
        # TODO: define the csv path and read the file
        # TODO: store it to system, so that you can access it later
        system.M_data = M_data
    else:
        pass

    # TODO: get interpolation 
    M_interp = get_interp(system.M_data, t)

    # TODO: apply the M, which is in the system base of 100 MVA
    system.GENROU.set()     # why not use `system.GENROU.alter()`?


def get_interp(M_data, t):
    """
    Interpolation for M machine data at time t.
    
    Parameters
    ----------
    M_data : pandas.DataFrame
        Two column dataframe, first for `t` and the second for `M`
    t : float
        time at which the interpolation is done
    """
    pass

    # # Define the file path to the CSV containing disturbance data
    # csv_file_path = 'disturbance_data.csv'

    # # Assuming you have a model instance and want to apply data to its 'x' attribute
    # selected_model = system.devices['your_model_name']  # Replace 'your_model_name' with the actual model name
    # selected_attribute = 'x'

    # try:
    #     # Read the CSV file into a DataFrame
    #     df = pd.read_csv(csv_file_path)

    #     # Check if the 'attribute' exists in the DataFrame
    #     if selected_attribute not in df.columns:
    #         raise ValueError(f"'{selected_attribute}' does not exist in the CSV file.")

    #     # Find the value in the CSV corresponding to the current simulation time
    #     current_row = df[df['time'] <= t].iloc[-1]

    #     # Extract the value from the CSV
    #     value = current_row[selected_attribute]

    #     # Apply the value to the model
    #     selected_model.set_attribute(selected_attribute, value)

    #     print(f"Time: {t}, {selected_attribute}: {value} applied to the model.")

    # except Exception as e:
    #     print(f"Error: {e}")

    # pass
