def PDF_Control_Dougherty(roll_error, roll_error_rate):
    # initializing Control variables as per Dougherty1968
    Kphi = 62.3676  # Control parameter in joules/rad
    tau = 25  # time constant of controller zero in sec

    Control_Torque = Kphi * (tau * roll_error_rate + roll_error)

    return Control_Torque
