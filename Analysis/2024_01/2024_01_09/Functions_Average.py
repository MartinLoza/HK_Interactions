
import cooler
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
import h5py
import matplotlib.pyplot as plt
import seaborn as sns
import warnings

def GetResolutionsMcool(mcool_path):
    """
    Get the available resolutions from an mcool file.
    
    Parameters:
    - mcool_path: str, path to the mcool file
    
    Returns:
    - resolutions: list of int, available resolutions in the mcool file
    """
    # Open the mcool file
    with h5py.File(mcool_path, 'r') as tmp_mcool:
        # Get the available resolutions
        resolutions = tmp_mcool['resolutions'].keys()
        # Change the resolutions to a vector
        resolutions = [int(res) for res in resolutions]
        # Order the resolutions
        resolutions.sort()
    
    return resolutions

def SameResolutionMcools(mcools_path):
    """
    Compare the resolutions of multiple mcool files.
    
    Parameters:
    - mcool_files: list of str, paths to the mcool files
    
    Returns:
    - resolutions_equal: bool, True if all resolutions are equal, False otherwise
    """
    resolutions_list = []
    
    for mcool_file in mcools_path:
        resolutions = GetResolutionsMcool(mcool_file)
        resolutions_list.append(resolutions)
    
    resolutions_equal = all(resolutions_list[i] == resolutions_list[i+1] for i in range(len(resolutions_list)-1))
    
    return resolutions_equal

def ScaleByGroup(df, group_column):
    """
    Scale the 'log10_balanced' values in a DataFrame by group.

    Parameters:
    - df (pandas.DataFrame): The input DataFrame.
    - group_column (str): The column name used for grouping.

    Returns:
    - scaled_matrix (pandas.DataFrame): The DataFrame with scaled 'log10_balanced' values.

    """
    # Obtain the unique values of the group column
    group_values = df[group_column].unique()

    # Create a list to store the scaled matrices
    scaled_matrices = []

    # Loop through each group value
    for value in group_values:
        # Filter the DataFrame by group value
        tmp_data = df[df[group_column] == value].copy()

        # Scale the log10_balanced values
        tmp_data['scaled'] = scaler.fit_transform(tmp_data['log10_balanced'].values.reshape(-1, 1))

        # Append the scaled matrix to the list
        scaled_matrices.append(tmp_data)

    # Concatenate the scaled matrices
    scaled_matrix = pd.concat(scaled_matrices)

    return scaled_matrix

def FilterIdNCells(df, n_cells=5, column='count', plot_histogram=True):
    """
    Filters the DataFrame based on the number of elements per ID and returns the filtered DataFrame.

    Parameters:
    - df: DataFrame
        The input DataFrame containing the data.
    - n_cells: int, optional (default=5)
        The number of cells required for an interaction ID to be included in the filtered DataFrame.
    - column: str, optional (default='count')
        The column in the DataFrame used for counting the number of elements per ID.
    - plot_histogram: bool, optional (default=True)
        Whether to plot a histogram of interaction IDs with less than n_cells elements.

    Returns:
    - df: DataFrame
        The filtered DataFrame containing only the interaction IDs that appear in n_cells cell types.
    """
    # count the number of elements per id
    tmp = df.groupby('id')[column].count()
    # ids of interaction appearing in n cell types
    ids = tmp[tmp == n_cells].index
    # filter the data
    df = df[df['id'].isin(ids)]
    
    if plot_histogram:
        # histogram of interaction in less than n cell types
        plt.hist(tmp[tmp < n_cells], bins=10)
        plt.title(f"Distribution of Interaction IDs with Less than {n_cells} Elements\nColumn: {column}")
        plt.xlabel("Number of Elements")
        plt.ylabel("Frequency")
        plt.show()
    
    return df

# Function to replace NaN values with a value outside the range of the data
def ReplaceNaNs(df, column='scaled', replace_with=-10):
    """
    Replaces NaN values in a specified column of a DataFrame with a specified value.
    
    Args:
        df (pandas.DataFrame): The DataFrame containing the column with NaN values.
        column (str): The name of the column to replace NaN values in. Default is 'scaled'.
        replace_with (float): The value to replace NaN values with. Default is -10.
        
    Returns:
        pandas.DataFrame: The DataFrame with NaN values replaced in the specified column.
    """
        
    #if the value is true, throw a warning and change the nan_value to something more appropriate
    if df[column].min() < replace_with:
        # Change the value based on the minimum value of the scaled data
        replace_with = round(df[column].min() - 1)  
        warnings.warn(f"The replace_with is not outside the range of the scaled data. Changing the value to: {nan_value}")
        
    # replace the NaN values with the nan_value
    df[column] = df[column].fillna(replace_with)
    
    return df

def PreprocessPixels(metadata, chromosomes, resolution, test=False):
    """
    Preprocesses the pixels of a matrix for a given cell type. Perform log10 transformation and scaling by chromosome.

    Args:
        metadata (DataFrame): The metadata containing information about the cell types and mcool file paths.
        chromosomes (list): The list of chromosomes to include in the analysis.
        resolution (int): The resolution of the mcool.

    Returns:
        DataFrame: The preprocessed and scaled matrix.

    """
    # List to store scaled matrices
    scaled_matrices = []

    if test:
        cell_index = [0]
    else:
        cell_index = list(range(len(metadata)))

    # Loop through each cell type, e.g. each row in the metadata
    for i in cell_index:
        # Get the current cell type
        cell_type = metadata['Biosource'][i]
        # Get the current mcool file path
        mcool_file = metadata['local_path'][i]
           
        # Load the mcool file
        cool_data = cooler.Cooler(mcool_file + "::resolutions/" + str(resolution))
            
        # Load the balanced matrix
        cool_matrix = cool_data.matrix(balance=True, as_pixels=True, join=True)[:]

        # Remove pixels of chromosomes not included in the main analysis
        cool_matrix = cool_matrix[cool_matrix['chrom1'].isin(chromosomes)].reset_index(drop=True)
        cool_matrix = cool_matrix[cool_matrix['chrom2'].isin(chromosomes)].reset_index(drop=True) 
            
        # remove inter-chromosomal interactions
        cool_matrix = cool_matrix[cool_matrix['chrom1'] == cool_matrix['chrom2']].reset_index(drop=True) 

        # Add a column with the cell type
        cool_matrix['cell_type'] = cell_type

        # Calculate the log10 of balanced counts
        cool_matrix['log10_balanced'] = cool_matrix.groupby('chrom1')['balanced'].transform(lambda x: np.log10(x))
            
        # Scale the log10 transformed values by chromosome using scale_by_group function
        scaled_matrix = ScaleByGroup(cool_matrix, 'chrom1')
            
        # Append the scaled matrix to the list
        scaled_matrices.append(scaled_matrix)

    # Concatenate the scaled matrices
    scaled_matrix = pd.concat(scaled_matrices)
        
    return scaled_matrix

def GetAverageMatrix(scaled_data, metadata,  plot_histogram=True, replaceNaN_with=-10, filter_nCells=6):

    # Let's add an interaction_id, then group by interaction (or loop by interaction) to average the scaled data.
    # This can be a bit time consuming for high resolution bins.
        
    # add interaction id
    scaled_data['id'] = scaled_data['chrom1'].astype(str) + '_' + scaled_data['start1'].astype(str) + '_' + scaled_data['chrom2'].astype(str) + '_' + scaled_data['start2'].astype(str)

    #Filter interactions that doesn't appear in all cell types using the local function
    scaled_data = FilterIdNCells(df = scaled_data, n_cells=filter_nCells, column='count', plot_histogram=plot_histogram)

    # Replace NaN values with a value outside the range of the data
    scaled_data = ReplaceNaNs(scaled_data, column='scaled', replace_with=replaceNaN_with)
    # short test to make sure everything is working
    test = scaled_data['scaled'].isna().sum()
    if test > 0:
        print(f"There are {test} NaN values in the scaled column after replacing them with an arbitrary value.")
        # stop the loop if there are still NaN values
        # break

    # Group by interaction id and average the scaled values
    mean_values = scaled_data.groupby('id')['scaled'].mean()
    # Convert the mean_values series to a dataframe and reset the index
    mean_values = mean_values.reset_index()
    # Rename the columns
    mean_values.columns = ['id', 'mean']

    # Create a new dataframe to store the mean values
    # we can filter the data for only one cell type and remove the unnecessary columns
    data_mean = scaled_data[scaled_data['cell_type'] == metadata['Biosource'][0]].copy()
    # select only the columns we need
    data_mean = data_mean[['id', 'chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2']]

    # Merge the mean_values with the data_mean dataframe
    data_mean = data_mean.merge(mean_values, on='id', how='left')

    # Select the columns we need
    data_mean = data_mean[['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'mean']]

    return data_mean
      