from torch.utils.data import Dataset, DataLoader, TensorDataset
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.preprocessing import LabelEncoder
from tqdm import tqdm


def get_data_loader(X, y, spatial):

    torch.set_num_threads(16)

    # Simulate dataset
    np.random.seed(0)

    n_samples, n_features  = X.shape

    # Create data loaders
    train_dataset = TensorDataset(torch.from_numpy(X).float(), torch.from_numpy(spatials).float(), torch.from_numpy(y))
    #test_dataset = TensorDataset(torch.from_numpy(X_test), torch.from_numpy(y_test))
    train_loader = DataLoader(train_dataset, batch_size=4096, shuffle=True)
    test_loader = DataLoader(train_dataset, batch_size=4096, shuffle=False)
    
    return train_dataset, train_loader
