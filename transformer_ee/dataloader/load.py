import numpy as np
import pandas as pd
import torch
from transformer_ee.dataloader.pd_dataset import Normalized_pandas_Dataset_with_cache

def get_sample_indices(sample_size: int, config) -> tuple:
    """
    A function to get the indices of the samples

    sample_size:    the number of samples
    config:         the configuration dictionary
    return:         the indices of train, validation and test sets
    """
    seed = config["seed"]
    _indices = np.arange(sample_size)
    np.random.seed(seed)
    np.random.shuffle(_indices)

    test_size = config["test_size"]
    valid_size = config["valid_size"]

    # test_size and valid_size can be either int or float
    if isinstance(test_size, float):
        test_size = int(sample_size * test_size)
    if isinstance(valid_size, float):
        valid_size = int(sample_size * valid_size)

    train_indices = _indices[: sample_size - valid_size - test_size]
    valid_indices = _indices[sample_size - valid_size - test_size : sample_size - test_size]
    test_indices = _indices[sample_size - test_size :]

    print("train indices size:\t", len(train_indices))
    print("valid indices size:\t", len(valid_indices))
    print("test  indices size:\t", len(test_indices))

    return train_indices, valid_indices, test_indices

def save_indices(indices, filename):
    """
    Save the indices to a text file.
    
    indices: The indices to save.
    filename: The file to save the indices in.
    """
    with open(filename, 'w') as f:
        for idx in indices:
            f.write(f"{idx}\n")

def get_train_valid_test_dataloader(config: dict):
    """
    A function to get the train, validation and test datasets
    Use the statistic of the training set to normalize the validation and test sets
    """
    df = pd.read_csv(config["data_path"])
    train_idx, valid_idx, test_idx = get_sample_indices(len(df), config)
    
    # Save the test indices
    save_indices(test_idx, config["save_path"]+'test_indices.txt')

    train_set = Normalized_pandas_Dataset_with_cache(
        config, df.iloc[train_idx].reset_index(drop=True, inplace=False)
    )
    valid_set = Normalized_pandas_Dataset_with_cache(
        config,
        df.iloc[valid_idx].reset_index(drop=True, inplace=False),
        weighter=train_set.weighter,
    )
    test_set = Normalized_pandas_Dataset_with_cache(
        config,
        df.iloc[test_idx].reset_index(drop=True, inplace=False),
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