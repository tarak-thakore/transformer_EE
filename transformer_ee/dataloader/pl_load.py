"""
Load the data
"""

import numpy as np
import polars as pl
import torch

from transformer_ee.dataloader.polars_dataset import (
    Normalized_Polars_Dataset_with_cache,
)
from .file_reader import get_polars_df_from_file


def get_sample_sizes(sample_size: int, config) -> tuple:
    """
    A function to get the indices of the samples

    sample_size:    the number of samples
    config:         the configuration dictionary
    return:         a tuple of three integers, the number of training samples, the number of validation samples, the number of test samples
    """

    test_size = config["test_size"]
    valid_size = config["valid_size"]

    # test_size and valid_size can be either int or float
    if isinstance(test_size, float):
        test_size = int(sample_size * test_size)
    if isinstance(valid_size, float):
        valid_size = int(sample_size * valid_size)

    print("train indicies size:\t", sample_size - test_size - valid_size)
    print("valid indicies size:\t", valid_size)
    print("test  indicies size:\t", test_size)

    return sample_size - test_size - valid_size, valid_size, test_size


def get_train_valid_test_dataloader(config: dict):
    """
    A function to get the train, validation and test datasets
    Use the statistic of the training set to normalize the validation and test sets
    """
    df = get_polars_df_from_file(config["data_path"])

    randomdf = df.sample(fraction=1.0, seed=config["seed"], shuffle=True)
    del df
    sizes = get_sample_sizes(randomdf.height, config)
    train_set = Normalized_Polars_Dataset_with_cache(
        config, randomdf.slice(offset=0, length=sizes[0])
    )
    valid_set = Normalized_Polars_Dataset_with_cache(
        config,
        randomdf.slice(offset=sizes[0], length=sizes[1]),
        weighter=train_set.weighter,
    )
    test_set = Normalized_Polars_Dataset_with_cache(
        config,
        randomdf.slice(offset=sizes[0] + sizes[1], length=sizes[2]),
        weighter=train_set.weighter,
    )

    train_set.statistic()

    train_set.normalize()
    valid_set.normalize(train_set.stat)
    test_set.normalize(train_set.stat)

    batch_size_train = config["batch_size_train"]
    batch_size_valid = config["batch_size_valid"]
    batch_size_test = config["batch_size_test"]

    trainloader = torch.utils.data.DataLoader(
        train_set,
        batch_size=batch_size_train,
        shuffle=True,
        num_workers=10,
    )

    validloader = torch.utils.data.DataLoader(
        valid_set,
        batch_size=batch_size_valid,
        shuffle=False,
        num_workers=10,
    )

    testloader = torch.utils.data.DataLoader(
        test_set,
        batch_size=batch_size_test,
        shuffle=False,
        num_workers=10,
    )

    return trainloader, validloader, testloader, train_set.stat
