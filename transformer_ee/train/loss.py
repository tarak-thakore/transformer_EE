"""
define loss function
"""

import torch

# shape of output and target: (batch_size, output_dim)
# shape of weight: (batch_size, 1)


# base loss functions
def MAE_loss(output, target, weight=None):
    """
    mean absolute error loss
    """
    if weight is None:
        return torch.mean(torch.abs(output - target))
    #print("output = ",output," and target = ", target,", so we have torch.mean(torch.abs(output - target)) = ",torch.mean(torch.abs(output - target)))
    return torch.mean(weight * torch.abs(output - target))


def MSE_loss(output, target, weight=None):
    """
    mean squared error loss
    """
    # Note: torch.mean() returns the mean value of all elements in the input tensor, which is a scalar value.
    if weight is None:
        return torch.mean((output - target) ** 2)
    return torch.mean(weight * (output - target) ** 2)

#This one is experimental...added by J. L. Barrow
def MCE_loss(output, target, weight=None):
    """
    mean cubic error loss
    """
    # Note: torch.mean() returns the mean value of all elements in the input tensor, which is a scalar value.
    if weight is None:
        return torch.mean((output - target) ** 3)
    return torch.mean(weight * (output - target) ** 3)


def MAPE_loss(output, target, weight=None):
    """
    mean absolute percentage error loss, similar to L1 loss
    """
    if weight is None:
        return torch.mean(torch.abs(output - target) / target)
    return torch.mean(weight * torch.abs(output - target) / target)


# add base loss functions to loss_function dictionary
loss_function = {
    "mean squared error": MSE_loss,
    "mean cubic error": MCE_loss,
    "mean absolute error": MAE_loss,
    "mean absolute percentage error": MAPE_loss,
}


# set up a base loss function by name
def get_loss_function(loss_function_name, output, target, weight=None):
    """
    get loss function
    """
    return loss_function[loss_function_name](output, target, weight)

######
# complex loss functions utilize base loss functions

def linear_combination_loss(output, target, weight=None, **kwargs):
    """
    linear combination of base loss functions
    coefficients, base_loss_names should have the same length, which is the number of output variables
    e.g. kwargs = {"coefficients": [0.5, 0.5], "base_loss_names": ["mean squared error", "mean absolute error"]}
    """
    if "base_loss_names" not in kwargs or "coefficients" not in kwargs:
        raise ValueError("base_loss_names and coefficients must be provided in kwargs")

    if len(kwargs["base_loss_names"]) != len(kwargs["coefficients"]):
        raise ValueError(
            "base_loss_names and coefficients must have the same length\n",
            "len(base_loss_names):",
            len(kwargs["base_loss_names"]),
            "\nlen(coefficients):",
            len(kwargs["coefficients"]),
        )

    base_loss_names = kwargs["base_loss_names"]
    coefficients = kwargs["coefficients"]
    linear_loss = 0
    for i in range(len(base_loss_names)):
        linear_loss += coefficients[i] * loss_function[base_loss_names[i]](
            output[:, i], target[:, i], torch.squeeze(weight)
        )
    return linear_loss

#Highly experimental--added by J. L. Barrow
def mass_constraint_loss(output, target):  #, weight=None, *kwargs):
    """
    Attempt to make a loss function which contrains the mass of the neutrino to zero; for use only with (E,px,py,pz)
    variables as targets
    """
    mass_neutrino_output = output[:,0] ** 2 - output[:,1] ** 2 - output[:,2] ** 2 - output[:,3] ** 2
    mass_neutrino_target = target[:,0] ** 2 - target[:,1] ** 2 - target[:,2] ** 2 - target[:,3] ** 2
    mass_loss = mass_neutrino_output - mass_neutrino_target
    return mass_loss